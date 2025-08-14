#include "fd.h"

int main(long argc, char *argv[])
{
  FILE *ppsi, *peval; par_st par;  long j, i, a, indexfirsthomo, thomo, tlumo,flags=0; 
  zomplex *potq, *psi, *phi;
  lng_st ist; fftw_plan_loc *planfw, *planbw; fftw_complex *fftwpsi;
  double *ksqr, *vx, *vy, *vz, *poth, *potl, *rx, *ry, *rz, tmp;
  double *psibe, *eval, *de, *Cbs, *Hbs, *Ebs;
  double *psiai, *evalai, *sigeai;
  double *vkijb, *vkcab;

  /*************************************************************************/
  writeCurrentTime(stdout);
  writeSeparation(stdout);
  init_size(argc,argv,&par,&ist);
  
  /*************************************************************************/  

  fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.ngrid*ist.nthreads);
  potq  = (zomplex*)calloc(ist.ngrid,sizeof(zomplex));
  poth = (double*)calloc(ist.ngrid*ist.nthreads,sizeof(double));
  potl = (double*)calloc(ist.ngrid,sizeof(double));
  ksqr = (double*)calloc(ist.ngrid,sizeof(double));
  vx = (double*)calloc(ist.nx,sizeof(double));
  vy = (double*)calloc(ist.ny,sizeof(double));
  vz = (double*)calloc(ist.nz,sizeof(double));
  rx = (double*)calloc(ist.natom,sizeof(double));
  ry = (double*)calloc(ist.natom,sizeof(double));
  rz = (double*)calloc(ist.natom,sizeof(double));
  psi  = (zomplex*)calloc(ist.ngrid,sizeof(zomplex));
  phi  = (zomplex*)calloc(ist.ngrid,sizeof(zomplex));
  Cbs = (double*)calloc(ist.msbs2*ist.numExcitons,sizeof(double));
  Hbs = (double*)calloc(ist.msbs2*ist.msbs2,sizeof(double));
  Ebs = (double*)calloc(ist.numExcitons,sizeof(double));
  
  /**************************************************************************/

  init(vx,vy,vz,ksqr,potl,rx,ry,rz,Cbs,Hbs,&par,&ist);

  /**************************************************************************/

  /*** initialization for the fast Fourier transform ***/
  fftw_plan_with_nthreads(ist.nthreads);
  planfw = (fftw_plan_loc *) calloc(ist.nthreads, sizeof(fftw_plan_loc));
  planbw = (fftw_plan_loc *) calloc(ist.nthreads, sizeof(fftw_plan_loc));
  for (i = 0; i < ist.nthreads; i++){ 
    planfw[i] = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,&fftwpsi[i*ist.ngrid],&fftwpsi[i*ist.ngrid],FFTW_FORWARD,flags);
    planbw[i] = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,&fftwpsi[i*ist.ngrid],&fftwpsi[i*ist.ngrid],FFTW_BACKWARD,flags);
  }

  /**************************************************************************/

  init_pot(vx, vy, vz, potq, par, ist, planfw[0], planbw[0], &fftwpsi[0]);

  /**************************************************************************/

  gauss_test(vx,vy,vz,potq,&poth[0],par,ist,planfw[0],planbw[0],&fftwpsi[0]);

  /**************************************************************************/

  get_energy_range(vx,vy,vz,ksqr,potl,&par,ist,planfw[0],planbw[0],&fftwpsi[0]);
  
  /**************************************************************************/

  eval = (double *) calloc(ist.homoIndex+1, sizeof(double)); 
  de = (double *) calloc(ist.homoIndex+1, sizeof(double)); 
  peval = fopen("eval.par" , "r");
  for (i = 0; i < ist.homoIndex+1; i++)
    fscanf (peval,"%ld %lg %lg",&a,&eval[i],&de[i]);
  fclose (peval);
  
  for (thomo = 0, i = ist.homoIndex; i >= 0; i--){
    if (de[i] < par.deps) thomo++;
    if (thomo == ist.numBandEdgeHoles) {
      indexfirsthomo = i;
      break;
    }
  }
  free(eval);
  free(de);

  psibe = (double *) calloc(ist.numBandEdgeStates*ist.ngrid, sizeof(double)); // band-edge, non-interacting states
  eval = (double *) calloc(ist.numBandEdgeStates, sizeof(double)); // band-edge, non-interacting energies 
  de = (double *) calloc(ist.numBandEdgeStates, sizeof(double)); // band-edge, non-interacting sigma values
  peval = fopen("eval.par" , "r");
  ppsi = fopen("psi.par" , "r");
  // skip over states not within boltzEnergyRange of the homo
  for (i = 0; i < indexfirsthomo; i++) {
    fscanf(peval, "%ld %lg %lg", &a, &tmp, &tmp);
    fread(&poth[0], sizeof(double), ist.ngrid, ppsi);
  }
  // store hole states within boltzEnergyRange of the homo 
  for (a = thomo = 0; (thomo < ist.numBandEdgeHoles) && (a <= ist.homoIndex); a++) {
    fscanf(peval, "%ld %lg %lg", &a, &eval[thomo], &de[thomo]);
    fread(&poth[0], sizeof(double), ist.ngrid,ppsi);
    if (de[thomo] < par.deps) {
      for (j = 0; j < ist.ngrid; j++) psibe[thomo*ist.ngrid+j] = poth[j];
      thomo++;
    }
  }
  // skip over states not that are not eigenstates in the gap
  for (i = ist.homoIndex; i < ist.lumoIndex-1; i++) {
    fscanf(peval, "%ld %lg %lg", &a, &tmp, &tmp);
    fread(&poth[0], sizeof(double), ist.ngrid, ppsi);
  }
  // store electron states within boltzEnergyRange of the lumo 
  for (tlumo = 0; tlumo < ist.numBandEdgeElectrons; ) {
    fscanf(peval, "%ld %lg %lg", &a, &eval[thomo+tlumo], &de[thomo+tlumo]);
    fread(&poth[0], sizeof(double), ist.ngrid, ppsi);
    if (de[thomo+tlumo] < par.deps) {
      for (j = 0; j < ist.ngrid; j++) psibe[(thomo+tlumo)*ist.ngrid+j] = poth[j];
      tlumo++;
    }
  }
  fclose(peval);
  fclose(ppsi);
  
  ist.homoIndex = ist.numBandEdgeHoles-1;
  ist.lumoIndex = ist.numBandEdgeHoles;

  ist.numBiexcitons = calc_bse_biexcitonic_states(Ebs, Cbs, Hbs, ist, par);
  printf("Number of initial biexcitonic states = %ld\n", ist.numBiexcitons);
  printf("Energy conservation window has DeltaE = %.6f\n", par.DeltaE); fflush(0);

  peval = fopen("eval.dat" , "w");
  for (i = 0; i < ist.numBandEdgeStates; i++) fprintf(peval, "%ld %.8f %.8f\n", i, eval[i], de[i]);
  fclose(peval);

  /**************************************************************************/
  
  if ((psiai = (double *) calloc(ist.ngrid*ist.ms, sizeof(double)))==NULL)nerror("psiai");
  printf("ist.ngrid*ist.ms = %ld\n", ist.ngrid*ist.ms); 
  evalai = (double *) calloc(ist.ms, sizeof(double));
  sigeai = (double *) calloc(ist.ms, sizeof(double));
  generate_filter_states(psiai,evalai,sigeai,eval,ksqr,potl,vx,vy,vz,&ist,&par,planfw,planbw,fftwpsi);
 
  /**************************************************************************/

  vkijb = (double *) calloc(ist.numBandEdgeHoles*ist.numBandEdgeHoles*ist.numBandEdgeElectrons*ist.itot, sizeof(double));
  vkcab = (double *) calloc(ist.numBandEdgeHoles*ist.numBandEdgeElectrons*ist.numBandEdgeElectrons*ist.atot, sizeof(double));
  generate_coulomb_matrix_elements(vkijb,vkcab,psiai,evalai,sigeai,psibe,potq,poth,ist,par,planfw,planbw,fftwpsi);
  
  /**************************************************************************/

  calculate_auger(Cbs,Ebs,vkijb,vkcab,eval,evalai,sigeai,ist,par);
    
  /**************************************************************************/
  free(potq); free(potl);  free(psibe); free(psiai); free(eval);
  free(ksqr); free(vx); free(vy);  free(vz); free(vkcab); free(vkijb);
  free(rx); free(ry); free(rz); free(poth); free(evalai); free(sigeai);
  free(Ebs); free(Cbs); free(Hbs);
  for (i = 0; i < ist.nthreads; i++){ 
    fftw_destroy_plan(planfw[i]);
    fftw_destroy_plan(planbw[i]);
  }
  fftw_free(fftwpsi);

  /**************************************************************************/
  writeSeparation(stdout);
  writeCurrentTime(stdout);

  exit(0);
}
