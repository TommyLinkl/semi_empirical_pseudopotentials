#include "fd.h"

/*****************************************************************************/

void init_size(long argc, char *argv[], par_st *par, lng_st *ist)
{
  FILE* pf = fopen("input.par" , "r");
  double evalloc, deloc, *eval, *de, excEnergy, tmp;
  long a;
  int i, ieof, nval;
 
  fscanf (pf,"%ld",&ist->nx);  /*** number of grid point in x ***/
  fscanf (pf,"%ld",&ist->ny);  /*** number of grid point in y ***/
  fscanf (pf,"%ld",&ist->nz);  /*** number of grid point in z ***/
  fscanf (pf,"%ld",&ist->nmc);  /*** number of iteration ***/ /*** ns in filter code, num of filter cycles - tommy***/
  fscanf (pf,"%ld",&ist->two);  /*** num of states per filter - tommy ***/
  fscanf (pf,"%ld",&ist->nc);    /*** length of newton interpolation ***/
  fscanf (pf,"%lg",&par->DeltaE); /*** DeltaE window used in enforcing energy conservation ***/
  fscanf (pf,"%ld",&ist->seed);    
  fscanf (pf,"%ld",&ist->nthreads);
  fscanf (pf,"%lg",&par->deps); /*** are you an eigenstate? ***/
  fscanf (pf,"%lg",&par->temp); /*** temperature for boltzmann averaging over initial states ***/
  fscanf (pf,"%s", ist->crystalStructure); /*** wurtzite or zincblende ***/
  fscanf (pf,"%s", ist->outmostMaterial); /*** CdS, CdSe, InP, InAs, alloyInGaP, alloyInGaAs ***/
  fclose(pf);

  pf = fopen("conf.par" , "r");
  fscanf(pf,"%ld",&ist->natom);
  fclose(pf);
  
  /*** set default parameters - no longer specified in input.par - John ***/
  ist->npot = 8192;	/*** size of pseudopotential files ***/
  par->Elmin = -1;	/*** minimum energy ***/
  par->Elmax = 1;	/*** maximum energy ***/
  par->Ekinmax = 10;	/*** kinetic energy cutoff ***/
  // ist->two = 20;   /*** num of states per filter is hard coded - tommy ***/
  ist->ms = ist->two * ist->nmc;	/*** number of states per filter ***/
  par->kbT = KB*par->temp/AUTOEV;
  par->boltzEnergyRange = 3.0*par->kbT;

  ist->nx_1 = 1.0 / (double)(ist->nx);
  ist->ny_1 = 1.0 / (double)(ist->ny);
  ist->nz_1 = 1.0 / (double)(ist->nz);
  ist->ngrid = ist->nx * ist->ny * ist->nz;
  ist->ngrid_1 = 1.0 / (double)(ist->ngrid);
  ist->natomtype = 15;
  
  printf("nx = %ld  ny = %ld  nz = %ld npot = %ld\n",
	  ist->nx, ist->ny, ist->nz, ist->npot);
  printf("ms = %ld  nmc = %ld natom = %ld nc = %ld nthreads = %ld\n",
	  ist->ms, ist->nmc, ist->natom, ist->nc, ist->nthreads);
  printf("Elmin = %g  Elmax = %g  Ekinmax = %g\n",
	  par->Elmin,par->Elmax,par->Ekinmax);
  printf("Temperature = %.3f kbT = %.8f boltzEnergyRange = %.8f\n", 
	  par->temp, par->kbT, par->boltzEnergyRange);
  printf("Crystal structure type is set to %s\n", ist->crystalStructure);
  printf("Outmost material type is set to %s\n", ist->outmostMaterial);

  ist->homoIndex = ist->lumoIndex = ist->totalHomo = ist->totalLumo = 0;
  pf = fopen("eval.par" , "r");
  for (i = ieof = 0; ieof != EOF; i++) {
    ieof = fscanf (pf,"%ld %lg %lg",&a,&evalloc,&deloc);
    if (deloc < par->deps && evalloc < -0.2) ist->homoIndex = i;
    if (deloc < par->deps && evalloc < -0.2) ist->totalHomo++; 
    if (deloc < par->deps && evalloc > -0.2 && ieof != EOF) ist->totalLumo++; 
  }
  fclose(pf);
  nval = i - 1; // the total number of states in original eval.par/psi.par
  pf = fopen("eval.par" , "r");
  for (i = 0; i <= ist->homoIndex; i++) fscanf (pf,"%ld %lg %lg",&a,&evalloc,&deloc);
  for (i = ist->homoIndex+1; i < nval; i++) {
    fscanf (pf,"%ld %lg %lg",&a,&evalloc,&deloc);
    if (deloc < par->deps) {
      ist->lumoIndex = i;
      break;
    }
  }
  fclose(pf);

  // Can just do this since BSE-based AR calculation requires BSE to be run first;
  // thus, the psi.par already contains only the band-edge eigenstates 
  // that were used in the BSE calculation
  ist->numBandEdgeHoles = ist->totalHomo;
  ist->numBandEdgeElectrons = ist->totalLumo;
  ist->numBandEdgeStates = ist->numBandEdgeHoles + ist->numBandEdgeElectrons;
  ist->msbs = ist->numBandEdgeHoles + ist->numBandEdgeElectrons;
  ist->msbs2 = ist->numBandEdgeHoles * ist->numBandEdgeElectrons;

  printf("Total number of hole eigenstates = %ld\n",ist->numBandEdgeHoles);
  printf("Total number of electron eigenstates = %ld\n",ist->numBandEdgeElectrons);
  printf("Total number of non-interacting states used in the BSE = %ld\n", ist->numBandEdgeStates); 
  printf("Number of hole eigenstates times number of electron eigenstates = %ld\n", ist->msbs2);
  printf("The index of the homo state = %ld\n",ist->homoIndex);
  printf("The index of the lumo state = %ld\n",ist->lumoIndex);

  eval = (double *)calloc(nval,sizeof(double)); 
  de = (double *)calloc(nval,sizeof(double)); 
  pf = fopen("eval.par" , "r");
  for (i = 0; i < nval; i++){
    fscanf (pf,"%ld %lg %lg", &a, &eval[i], &de[i]);
  }
  fclose (pf);

  par->homoEnergy = eval[ist->homoIndex];
  par->lumoEnergy = eval[ist->lumoIndex];
  par->Egap = par->lumoEnergy - par->homoEnergy;
  printf("Homo energy = %.10f\n", par->homoEnergy);
  printf("Lumo energy = %.10f\n", par->lumoEnergy);
  printf("Energy gap  = %.10f %.8f\n", par->Egap, par->Egap*AUTOEV);

  free(eval);  free(de);

  // Read in the minimum initial biexcitonic state (product of two lowest excitonic states) 
  // energy (stored in par->minInitE) and calculate the  number of excitonic states 
  // within boltzEnergyRange (= 3kbT) of the lowest energy excitonic state (ist->numExcitons) 
  ist->numExcitons = 0;
  pf = fopen("exciton.par" , "r");
  for (i = ieof = 0; ieof != EOF; i++) {
    ieof = fscanf (pf,"%ld %lg %lg %lg %lg", &a, &excEnergy, &tmp, &tmp, &tmp);
    if (i == 0) { 
      par->Ex = excEnergy; // energy of lowest excitonic state
      ist->numExcitons++;
    }
    // make sure excitonic state is in desired energy range to be included as a possible 
    // initial state (assumed 1 of the initial excitons is the lowest energy exciton)
    else if (excEnergy <= (par->Ex + par->boltzEnergyRange)) ist->numExcitons++;
    else break; // all other excitonic states are of too high energy
  }
  fclose(pf);

  par->minInitE = 2.0*par->Ex;
  par->maxInitE = 2.0*par->Ex + par->boltzEnergyRange;
  printf("Ground state exciton has energy = %.10f %.8f\n", par->Ex, par->Ex*AUTOEV);
  printf("Number of excitons within boltzEnergyRange of the ground state energy = %ld\n", ist->numExcitons);
  printf("Min energy of the initial ground state biexciton =  %.10f\n", par->minInitE);
  printf("Max allowed energy of the initial biexciton =  %.10f\n", par->maxInitE);

  return;
}

/*****************************************************************************/

void init(double *vx,double *vy,double *vz,double *ksqr,double *potl,double *rx,double *ry,double *rz,double *Cbs,double *Hbs,par_st *par,lng_st *ist)
{
  FILE *pf, *pfs; 
  atm_st *atm;
  long jexc, jx, jy, jz, jyz, jxyz, iatom, ntot, jp, *npot, nn, flags=0, ii;
  double del, mx, my, mz, xd, yd, zd, dx, dy, dz, *ksqrx, *ksqry, *ksqrz;
  double *vr, *potatom, *dr, sum, rex, potEx, *a4Params, *a5Params, *volRef, *strainScale, *tetrahedronVol; 
  pot_st ppar;
  vector *atomNeighborList;
  int crystalStructureInt, outmostMaterialInt;

  if ((ksqrx  = (double*)calloc(ist->nx,sizeof(double)))==NULL)nerror("ksqrx");
  if ((ksqry  = (double*)calloc(ist->ny,sizeof(double)))==NULL)nerror("ksqry");
  if ((ksqrz  = (double*)calloc(ist->nz,sizeof(double)))==NULL)nerror("ksqrz");
  if ((atm = (atm_st*)calloc(ist->natom,sizeof(atm_st)))==NULL)nerror("atm");
 
// Read in the external potential and store it in ppar if it exists
  if ( access( "expot.par", F_OK) != -1 ) {
  pf = fopen("expot.par", "r");
  fscanf (pf,"%lg",&ppar.e);  /*** external potential parameter: Energy Scale ***/
  fscanf (pf,"%lg",&ppar.x0); /*** external potential parameter: Ef ***/
  fscanf (pf,"%lg",&ppar.t);  /*** external potential parameter: kT ***/
  fclose(pf);
  printf("expot.par was detected\n");
  printf("External potential energy scale, e = % .6f\n", ppar.e);
  printf("External potential Fermi radius, x0 = % .6f\n", ppar.x0);
  printf("External potential steepness (Temp), kT = % .6f\n", ppar.t);
  } 
  else {
  ppar.e = ppar.x0 = 0.0; // setting ppar.e = 0 results in ext pot being 0
  ppar.t = 1.0;
  printf("No expot.par file detected in cwd - setting external potential to 0!\n");
  }
  
  /*** read the pasivated nanocrystal configuration ***/
  pf = fopen("conf.par" , "r");
  fscanf(pf,"%ld",&ntot);
  assert(fabs((double)(ntot - ist->natom)) < 1.0e-15);
  read_conf(rx,ry,rz,atm,ntot,pf);
  fclose(pf);
  
  xd = rint(0.5 * get_dot_ligand_size_z(rx,ntot) + 5.0);
  yd = rint(0.5 * get_dot_ligand_size_z(ry,ntot) + 5.0);
  zd = rint(0.5 * get_dot_ligand_size_z(rz,ntot) + 5.0);
  printf("Box (quadrant) dimensions: xd = %.2f yd = %.2f zd = %.2f\n", xd, yd, zd);

  /***initial parameters for the pot reduce mass, etc. in the x direction ***/
  par->xmin = -xd;
  par->xmax = xd;
  mx = 1.0;
  par->dx  = (par->xmax - par->xmin) / (double)(ist->nx);
  par->dkx = TWOPI / ((double)ist->nx * par->dx);
  
  /***initial parameters for the pot reduce mass, etc. in the y direction ***/
  par->ymin = -yd;
  par->ymax = yd;
  my = 1.0;
  par->dy  = (par->ymax - par->ymin) / (double)(ist->ny);
  par->dky = TWOPI / ((double)ist->ny * par->dy);

  /***initial parameters for the pot reduce mass, etc. in the z direction ***/
  par->zmin = -zd;
  par->zmax = zd;
  mz = 1.0;
  par->dz  = (par->zmax - par->zmin) / (double)(ist->nz);
  par->dkz = TWOPI / ((double)ist->nz * par->dz);

  par->Vmin = 1.0e10;
  par->Vmax = -1.0e10;

  par->gamma = 7.0 / (2.0 * xd);
  par->gamma2 = sqr(par->gamma);

  par->dv = par->dx * par->dy * par->dz;
  par->dr = sqrt(sqr(par->dx) + sqr(par->dy) + sqr(par->dz));
  printf("Grid point spacing: dx = %.4f dy = %.4f dz = %.4f dv = %.6f dr = %.6f\n",
    par->dx, par->dy, par->dz, par->dv, par->dr);

  /***initializing the ksqr vectors ***/
  for (ksqrx[0] = 0.0, jx = 1; jx <= ist->nx / 2; jx++)
  ksqrx[jx] = (ksqrx[ist->nx-jx] = 0.5 * sqr((double)(jx) * par->dkx) *
    ist->nx_1 * ist->ny_1 * ist->nz_1 / mx);

  for (ksqry[0] = 0.0, jy = 1; jy <= ist->ny / 2; jy++)
  ksqry[jy] = (ksqry[ist->ny-jy] = 0.5 * sqr((double)(jy) * par->dky) *
    ist->ny_1 * ist->nx_1 * ist->nz_1 / my);

  for (ksqrz[0] = 0.0, jz = 1; jz <= ist->nz / 2; jz++)
  ksqrz[jz] = (ksqrz[ist->nz-jz] = 0.5 * sqr((double)(jz) * par->dkz) *
    ist->nz_1 * ist->nx_1 * ist->ny_1 / mz);

  par->Ekinmax *= (ist->ny_1 * ist->nx_1 * ist->nz_1);
  for (jz = 0; jz < ist->nz; jz++) for (jy = 0; jy < ist->ny; jy++){
  for (jyz = ist->nx * (ist->ny * jz + jy), jx = 0; jx < ist->nx; jx++){
    jxyz = jyz + jx;
    ksqr[jxyz] = ksqrx[jx] + ksqry[jy] + ksqrz[jz];
    if (ksqr[jxyz] > par->Ekinmax) ksqr[jxyz] = par->Ekinmax;
  }
  }
  free(ksqrx); free(ksqry);  free(ksqrz);

  /*** read pseudopotentials ***/
  fprintf(stdout, "Number of atom types = %ld\n", ist->natomtype); fflush(stdout);

  dr  = (double *) calloc(ist->natomtype, sizeof(double));
  vr  = (double *) calloc(ist->npot*ist->natomtype, sizeof(double));
  potatom = (double *) calloc(ist->npot*ist->natomtype, sizeof(double));
  npot = (long *) calloc(ist->natomtype, sizeof(long));
  a4Params = (double *) calloc(ist->natomtype, sizeof(double));
  a5Params = (double *) calloc(ist->natomtype, sizeof(double));
  read_pot(vr, potatom, npot, dr, atm, ist->npot, ist->natomtype, a4Params, a5Params);

  /*** calculate strainScale for atomic pseudopotential ***/
  if ((atomNeighborList = (vector *) calloc(4*ntot, sizeof(vector))) == NULL) nerror("atomNeighborList");
  if ((strainScale = (double *) calloc(ntot, sizeof(double))) == NULL) nerror("strainScale");
  if ((volRef = (double *) calloc(ntot, sizeof(double))) == NULL) nerror("volRef");

  /*** Assign crystalStructureInt ***/
  if (! strcmp(ist->crystalStructure, "wurtzite")) {
    crystalStructureInt = 0;
  }
  else if (! strcmp(ist->crystalStructure, "zincblende")) {
    crystalStructureInt = 1;
  }
  else {
  printf("\n\nCrystal structure type %s not recognized -- the program is exiting!!!\n\n", ist->crystalStructure);
  fflush(stdout);
  exit(EXIT_FAILURE);
  }

  /*** Assign outmostMaterialInt ***/
  if (! strcmp(ist->outmostMaterial, "CdS")) {
    outmostMaterialInt = 0;
  }
  else if (! strcmp(ist->outmostMaterial, "CdSe")) {
    outmostMaterialInt = 1;
  }
  else if (! strcmp(ist->outmostMaterial, "InP")) {
    outmostMaterialInt = 2;
  }
  else if (! strcmp(ist->outmostMaterial, "InAs")) {
    outmostMaterialInt = 3;
  }
  else if (! strcmp(ist->outmostMaterial, "alloyInGaP")) {  //cation terminated surface only
    outmostMaterialInt = 4;
  }
  else if (! strcmp(ist->outmostMaterial, "alloyInGaAs")) {  // cation terminated surface only
    outmostMaterialInt = 5;
  }
  else {
    printf("\n\nOutmostMaterial type %s not recognized -- the program is exiting!!!\n\n", ist->outmostMaterial);
    fflush(stdout);
    exit(EXIT_FAILURE);
  }


  readNearestNeighbors(ntot, crystalStructureInt, atomNeighborList, volRef, outmostMaterialInt);

  calculateStrainScale(ntot, volRef, atm, atomNeighborList, a4Params, a5Params, strainScale);

  free(atomNeighborList); free(volRef); free(a4Params); free(a5Params);

  /***initializing the potential vector  ***/
  for (jx = 0, dx = par->xmin; jx < ist->nx; jx++, dx += par->dx) vx[jx] = dx;
  for (jy = 0, dy = par->ymin; jy < ist->ny; jy++, dy += par->dy) vy[jy] = dy;
  for (jz = 0, dz = par->zmin; jz < ist->nz; jz++, dz += par->dz) vz[jz] = dz;

  omp_set_dynamic(0);
  omp_set_num_threads(ist->nthreads);
#pragma omp parallel for private(dx,dy,dz,del,jy,jx,jyz,jxyz,sum,iatom)
  for (jz = 0; jz < ist->nz; jz++) {
  for (jy = 0; jy < ist->ny; jy++) {
    jyz = ist->nx * (ist->ny * jz + jy);
    for (jx = 0; jx < ist->nx; jx++) {
    rex = sqrt(vx[jx]*vx[jx]+vy[jy]*vy[jy]+vz[jz]*vz[jz]);
    potEx = expot(rex,ppar);
    jxyz = jyz + jx;
    for (sum = 0.0, iatom = 0; iatom < ntot; iatom++) {
       dx = vx[jx] - rx[iatom];
       dy = vy[jy] - ry[iatom];
       dz = vz[jz] - rz[iatom];
       del = sqrt(dx * dx + dy * dy + dz * dz);
       sum += interpolate(del,dr[atm[iatom].natyp],vr,potatom,ist->npot,npot[atm[iatom].natyp],atm[iatom].natyp, strainScale[iatom]);
    } 
    potl[jxyz] = sum+potEx;
    }
    // printf("potl[jxyz] = %.4f\n", potl[jxyz]); 
  }
  }

  // Print out the external potential as a function of the distance from the origin 
  pfs = fopen("pS.dat", "w");
  for (jz = 0; jz < ist->nz; jz++) {
  potEx=expot(fabs(vz[jz]),ppar);
  fprintf(pfs,"% .6f % .8f\n", vz[jz], potEx);
  }
  fclose(pfs);

  writeCubeFile(potl, par->xmin, par->ymin, par->zmin, par->dx, par->dy, par->dz, ist->nx, ist->ny, ist->nz);

  for (jxyz= 0; jxyz < ist->ngrid; jxyz++){
  if (par->Vmax < potl[jxyz]) par->Vmax = potl[jxyz];
  if (par->Vmin > potl[jxyz]) par->Vmin = potl[jxyz];
  }
  /*write_pot(vx,vy,vz,potl);*/
  printf("Potential energy statistics: dV = %.5f Vmin = % .5f Vmax = %.5f\n",
    par->Vmax-par->Vmin,par->Vmin,par->Vmax);

  par->dE = 0.5 * sqr(PIE) / (mx*par->dx*par->dx) +
  0.5 * sqr(PIE) / (my*par->dy*par->dy) +    
  0.5 * sqr(PIE) / (mz*par->dz*par->dz);
  printf("dT = %g\n", par->dE);
  
  /*** initialization for the fast Fourier transform ***/
  //(*planfw) = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_FORWARD,flags);
  //(*planbw) = fftw_plan_dft_3d(ist->nz,ist->ny,ist->nx,fftwpsi,fftwpsi,FFTW_BACKWARD,flags);
  //(*planfw) = fftw3d_create_plan(ist->nz,ist->ny,ist->nx,FFTW_FORWARD,flags);
  //(*planbw) = fftw3d_create_plan(ist->nz,ist->ny,ist->nx,FFTW_BACKWARD,flags);


  /*** read the BSE coefficients ***/
  pf = fopen ("BSEcoeff.par" , "r");
  for (jexc = 0; jexc < ist->numExcitons; jexc++) {
    for (jx = 0; jx < ist->msbs2; jx++) fscanf(pf,"%lg", &Cbs[jx+jexc*ist->msbs2]);
  }
  fclose(pf);

  pf = fopen ("HBSEmat.par" , "r");
  for (jx = 0; jx < ist->msbs2*ist->msbs2; jx++) fscanf(pf, "%lg", &Hbs[jx]);
  fclose(pf);

  free(dr); free(vr); free (atm);  free(npot); 
  free(potatom); free(strainScale);

  return;
}

/****************************************************************************/

/****************************************************************************/
/* read all neighbors of each atom */
void readNearestNeighbors(long nAtoms, int crystalStructure, vector *atomNeighbors, double *tetrahedronVolRef, int outmostMaterial) {
  
  FILE *pf;
  long iAtom, at_natyp, *neighbors_natyp; 
  char at[4], n1[4], n2[4], n3[4], n4[4];
  char strerror[200];
  int tmpi, numSe, numS, numCd, numAs, numIn, numP, numGa, iNeighbor;
  double tmpx, tmpy, tmpz;
  double *bondLengths; 
  if ((neighbors_natyp = (long*)calloc(4,sizeof(double)))==NULL)nerror("neighbors_natyp");
  if ((bondLengths = (double*)calloc(4,sizeof(double)))==NULL)nerror("bondLengths");
  
  if ( access("allNeighborBonds.par", F_OK) != -1 ) {
    
    pf = fopen("allNeighborBonds.par", "r");
    
    for (iAtom = 0; iAtom < nAtoms; iAtom++) {
      fscanf(pf, "%s ", at);
      if ((! strcmp(at, "P1")) || (! strcmp(at, "P2")) || (! strcmp(at, "P3")) || (! strcmp(at, "P4"))) {
  fscanf(pf, "%i %lg %lg %lg %s %i %lg %lg %lg",
         &tmpi, &tmpx, &tmpy, &tmpz,
         n1, &tmpi, &atomNeighbors[4*iAtom].x, &atomNeighbors[4*iAtom].y, &atomNeighbors[4*iAtom].z);
  tetrahedronVolRef[iAtom] = 0.;
      }

      else {
  fscanf(pf, "%i %lg %lg %lg %s %i %lg %lg %lg %s %i %lg %lg %lg %s %i %lg %lg %lg %s %i %lg %lg %lg",
         &tmpi, &tmpx, &tmpy, &tmpz,
         n1, &tmpi, &atomNeighbors[4*iAtom].x, &atomNeighbors[4*iAtom].y, &atomNeighbors[4*iAtom].z,
         n2, &tmpi, &atomNeighbors[4*iAtom+1].x, &atomNeighbors[4*iAtom+1].y, &atomNeighbors[4*iAtom+1].z,
         n3, &tmpi, &atomNeighbors[4*iAtom+2].x, &atomNeighbors[4*iAtom+2].y, &atomNeighbors[4*iAtom+2].z,
         n4, &tmpi, &atomNeighbors[4*iAtom+3].x, &atomNeighbors[4*iAtom+3].y, &atomNeighbors[4*iAtom+3].z);
  
  // adjust positions of passivation ligands
  if (! strcmp(n3, "P1")) {
    strcpy(n3, n2);
    atomNeighbors[4*iAtom+2] = retScaledVector(atomNeighbors[4*iAtom+2], 1.0/0.55);
  }
  if (! strcmp(n4, "P1")) {
    strcpy(n4, n2);
    atomNeighbors[4*iAtom+3] = retScaledVector(atomNeighbors[4*iAtom+3], 1.0/0.55);
  }
  if (! strcmp(n3, "P2") && ! strcmp(n4, "P2")) {
    atomNeighbors[4*iAtom+2] = retScaledVector(atomNeighbors[4*iAtom+2], 1./0.25);
    atomNeighbors[4*iAtom+3] = retScaledVector(atomNeighbors[4*iAtom+3], 1./0.25);
    // Below is just for Eran's old geometry
    // atomNeighbors[4*iAtom+2] = retScaledVector(atomNeighbors[4*iAtom+2], 1./0.30);
    // atomNeighbors[4*iAtom+3] = retScaledVector(atomNeighbors[4*iAtom+3], 1./0.30);
  }
  else if (! strcmp(n4, "P2")) {
    atomNeighbors[4*iAtom+3] = retScaledVector(atomNeighbors[4*iAtom+3], 1.0/0.30);
    // Below is just for Eran's old geometry
    // atomNeighbors[4*iAtom+3] = retScaledVector(atomNeighbors[4*iAtom+3], 1.0/0.40);
  }
  
  at_natyp = assign_atom_number(at);
  neighbors_natyp[0] = assign_atom_number(n1);
  neighbors_natyp[1] = assign_atom_number(n2);
  neighbors_natyp[2] = assign_atom_number(n3);
  neighbors_natyp[3] = assign_atom_number(n4);
  
  for (iNeighbor=0; iNeighbor<4; iNeighbor++) {
    bondLengths[iNeighbor] = 0.;
    // If neighbor is a passivation ligand, replace it with the corresponding semiconductor atom
    if ((neighbors_natyp[iNeighbor]==8) || (neighbors_natyp[iNeighbor]==9) || (neighbors_natyp[iNeighbor]==10) || (neighbors_natyp[iNeighbor]==11)) { 
      if ((outmostMaterial==0) && (at_natyp==0)) neighbors_natyp[iNeighbor]=7; // CdS, Center-Cd, Replace with S. 
      else if ((outmostMaterial==0) && (at_natyp==7)) neighbors_natyp[iNeighbor]=0; // CdS, Center-S, Replace with Cd. 
      else if ((outmostMaterial==0) && (at_natyp!=0) && (at_natyp!=7)) {
        sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp); 
        nerror(strerror);
      } 
      else if ((outmostMaterial==1) && (at_natyp==0)) neighbors_natyp[iNeighbor]=1; // CdSe, Center-Cd, Replace with Se. 
      else if ((outmostMaterial==1) && (at_natyp==1)) neighbors_natyp[iNeighbor]=0; // CdSe, Center-Se, Replace with Cd. 
      else if ((outmostMaterial==1) && (at_natyp!=0) && (at_natyp!=1)) {
        sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp);
        nerror (strerror);
      }
      else if ((outmostMaterial==2) && (at_natyp==2)) neighbors_natyp[iNeighbor]=16; // InP, Center-In, Replace with P. 
      else if ((outmostMaterial==2) && (at_natyp==16)) neighbors_natyp[iNeighbor]=2; // InP, Center-P, Replace with In. 
      else if ((outmostMaterial==2) && (at_natyp!=2) && (at_natyp!=16)) {
        sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp);
        nerror (strerror);
      }
      else if ((outmostMaterial==3) && (at_natyp==2)) neighbors_natyp[iNeighbor]=3; // InAs, Center-In, Replace with As. 
      else if ((outmostMaterial==3) && (at_natyp==3)) neighbors_natyp[iNeighbor]=2; // InAs, Center-As, Replace with In. 
      else if ((outmostMaterial==3) && (at_natyp!=2) && (at_natyp!=3)) {
        sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp);
        nerror (strerror);
      }
      else if ((outmostMaterial==4) && ((at_natyp==2)||(at_natyp==15))) neighbors_natyp[iNeighbor]=16; // Outmost: alloyInGaP; Center: In or Ga. Replace with P. 
      else if ((outmostMaterial==4) && (at_natyp==16)) neighbors_natyp[iNeighbor]=2; // Outmost: alloyInGaP; Center: P. Replace with In. This is a completely random choice. 
      else if ((outmostMaterial==4) && (at_natyp!=2) && (at_natyp!=15) && (at_natyp!=16)) {
        sprintf(strerror,"Outmost layer is input as %d, but atom type %ld is bonded to passivation ligands\n", outmostMaterial, at_natyp);
        nerror (strerror);
      }
      else if ((outmostMaterial==5) && ((at_natyp==2)||(at_natyp==15))) neighbors_natyp[iNeighbor]=3; // Outmost: alloyInGaAs; Center: In or Ga. Replace with As. 
      else if ((outmostMaterial==5) && (at_natyp==3)) neighbors_natyp[iNeighbor]=2; // Outmost: alloyInGaAs; Center: As. Replace with In. This is a completely random choice. 
      else if ((outmostMaterial==5) && (at_natyp!=2) && (at_natyp!=15)) {
        sprintf(strerror,"Outmost layer is input as %d, but the surface is not cation terminated.\n", outmostMaterial);
        nerror (strerror);
      }
    }
    bondLengths[iNeighbor] = retIdealBondLength(at_natyp, neighbors_natyp[iNeighbor], crystalStructure); 
    // printf("at_natyp = %d, iNeighbor = %d, neighbors_natyp[iNeighbor] = %d, bondLengths[iNeighbor]=%g \n", at_natyp, iNeighbor, neighbors_natyp[iNeighbor], bondLengths[iNeighbor]); 
  }
  tetrahedronVolRef[iAtom] = retRegularTetrahedronVolume(bondLengths[0], bondLengths[1], bondLengths[2], bondLengths[3]);
  // printf("The reference tetrahedron volume is %.7f \n", tetrahedronVolRef[iAtom]);  
      }
    }
    fclose(pf);
  }
  else {
    printf("\n\nNo allNeighborBonds.par file in current working directory -- the program is exiting!!!\n\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
  }

  free(neighbors_natyp);
  free(bondLengths);

  return;
}

/****************************************************************************/
/* compute strain scaling factor for each atom
 * strainScale = 1. + a4 * (Omega/Omega0 - 1),
 * where Omega is volume of tetrahedron formed by atom's nearest neighbors */

void calculateStrainScale(long nAtoms, double *tetrahedronVolRef, atm_st *atm, vector *atomNeighbors, double *a4Params, double *a5Params, double *strainScale) {

  FILE *pf;
  long iAtom;
  vector v1, v2, v3, v4;
  char at[4];
  double tmp2, tmp3, *tetrahedronVol, strain;
  int tmp1;

  if ((tetrahedronVol = (double *) calloc(nAtoms, sizeof(double))) == NULL) nerror("tetrahedronVol");

  if ( access("strain.par", F_OK) != -1 ) {
  
    printf("Reading strain scale from strain.par...\n");    
    pf = fopen("strain.par", "r");
    for (iAtom = 0; iAtom < nAtoms; iAtom++) {
      fscanf(pf, "%i %s %lg %lg %lg %lg %lg", &tmp1, at, &tmp2, &tmp3, 
          &tetrahedronVol[iAtom], &tetrahedronVolRef[iAtom], &strainScale[iAtom]);
    }
    fclose(pf);

  }
  else {
  
    // TODO: parallelize
    for (iAtom = 0; iAtom < nAtoms; iAtom++) {

      if ((atm[iAtom].natyp == 8) || (atm[iAtom].natyp == 9)) { 
        strainScale[iAtom] = 1.;
      }
      else {

        v1 = retSubtractedVectors(atomNeighbors[4*iAtom], atomNeighbors[4*iAtom+3]);
        v2 = retSubtractedVectors(atomNeighbors[4*iAtom+1], atomNeighbors[4*iAtom+3]);
        v3 = retSubtractedVectors(atomNeighbors[4*iAtom+2], atomNeighbors[4*iAtom+3]);

        tetrahedronVol[iAtom] = fabs(retDotProduct(v1, retCrossProduct(v2, v3)))/6.;
        strain = (tetrahedronVol[iAtom]/tetrahedronVolRef[iAtom] - 1.);
        strainScale[iAtom] = (1. + a4Params[atm[iAtom].natyp]*strain + a5Params[atm[iAtom].natyp]*cube(strain));
      }

    }
  }

  pf = fopen("strain.dat","w");
  for (iAtom = 0; iAtom < nAtoms; iAtom++) {
    fprintf(pf, "%ld %s %.8f %.8f %.8f %.8f %.8f\n", iAtom, atm[iAtom].atyp, a4Params[atm[iAtom].natyp], a5Params[atm[iAtom].natyp], 
        tetrahedronVol[iAtom], tetrahedronVolRef[iAtom], strainScale[iAtom]);
  }
  fclose(pf);

  free(tetrahedronVol);

  return;
}

/****************************************************************************/

double retRegularTetrahedronVolume(double bondLength1, double bondLength2, double bondLength3, double bondLength4) {
  
  vector v1, v2, v3, v4; 

  v1.x=sqrt(8./9.);   v1.y=0.;   v1.z=-1./3.;
  v2.x=-sqrt(2./9.);   v2.y=sqrt(2./3.);   v2.z=-1./3.;
  v3.x=-sqrt(2./9.);   v3.y=-sqrt(2./3.);   v3.z=-1./3.;
  v4.x=0.;   v4.y=0.;   v4.z=1.;
  
  v1 = retScaledVector(v1, bondLength1);
  v2 = retScaledVector(v2, bondLength2);
  v3 = retScaledVector(v3, bondLength3);
  v4 = retScaledVector(v4, bondLength4);

  vector side1 = retSubtractedVectors(v1, v4); 
  vector side2 = retSubtractedVectors(v2, v4); 
  vector side3 = retSubtractedVectors(v3, v4); 

  double regularTetrahedronVolume = fabs(retDotProduct(side1, retCrossProduct(side2, side3)))/6.;
  return regularTetrahedronVolume; 
}

/****************************************************************************/

void init_pot(double *vx,double *vy,double *vz,zomplex *potq,par_st par,lng_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  FILE *pf;
  long jx, jy, jz, jyz, jxyz, sx, sy, sz;
  double dr, x2, y2, z2, *kx2, *ky2, *kz2, alpha, cosa, sina;
  double boxl = (double)(ist.nx) * par.dx;
  zomplex *potr, tmp;

  if ((kx2  = (double*)calloc(ist.nx,sizeof(double)))==NULL)nerror("kx2");
  if ((ky2  = (double*)calloc(ist.ny,sizeof(double)))==NULL)nerror("ky2");
  if ((kz2  = (double*)calloc(ist.nz,sizeof(double)))==NULL)nerror("kz2");
  if ((potr = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("potr");
  
  for (jx = 0; jx < ist.ngrid; jx++) potr[jx].re = potr[jx].im = 0.0;
  for (jz = 0; jz < ist.nz; jz++) {
    z2 = sqr(vz[jz]);
    for (jy = 0; jy < ist.ny; jy++) {
      y2 = sqr(vy[jy]);
      jyz = ist.nx * (ist.ny * jz + jy);
      for (jx = 0; jx < ist.nx; jx++) {
	    x2 = sqr(vx[jx]);
	    jxyz = jyz + jx;
	    dr = sqrt(x2 + y2 + z2);
	    if (dr < boxl) potr[jxyz].re = screenedcoulomb(dr,par.gamma);
      }
    }
  }
  
  memcpy(&fftwpsi[0],&potr[0],ist.ngrid*sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  memcpy(&potq[0],&fftwpsi[0],ist.ngrid*sizeof(potq[0]));
  
  for (kx2[0] = 0.0, jx = 1; jx <= ist.nx / 2; jx++)
    kx2[jx] = (kx2[ist.nx-jx] = sqr((double)(jx) * par.dkx));
  for (ky2[0] = 0.0, jy = 1; jy <= ist.ny / 2; jy++)
    ky2[jy] = (ky2[ist.ny-jy] = sqr((double)(jy) * par.dky));
  for (kz2[0] = 0.0, jz = 1; jz <= ist.nz / 2; jz++)
    kz2[jz] = (kz2[ist.nz-jz] = sqr((double)(jz) * par.dkz));

  pf = fopen("qSpaceCoulombPot.dat", "w");
  for (sz = 1.0, jz = 0; jz < ist.nz; jz++, sz = -sz) {
    z2 = kz2[jz];
    for (sy = 1.0, jy = 0; jy < ist.ny; jy++, sy = -sy) {
      y2 = ky2[jy];
      jyz = ist.nx * (ist.ny * jz + jy);
      for (sx = 1.0, jx = 0; jx < ist.nx; jx++, sx = -sx) {
      	x2 = kx2[jx];
      	jxyz = jyz + jx;
      	alpha = PIE * (double)(jx + jy + jz + ist.ngrid / 2);
      	cosa = cos(alpha);
      	sina = sin(alpha);

      	tmp.re = potq[jxyz].re;
      	tmp.im = potq[jxyz].im;

      	potq[jxyz].re = (tmp.re * cosa - tmp.im * sina) * par.dv;
      	potq[jxyz].im = (tmp.re * sina + tmp.im * cosa) * par.dv;
      	
        double  q2 = x2 + y2 + z2 + 1.0e-20;
        /*potq[jxyz].re += (FOURPI / (q2 + par.gamma2));*/
        potq[jxyz].re += (FOURPI * (1.0 - exp(-0.25* q2 / par.gamma2)) / q2);


      	potq[jxyz].re *= ist.ngrid_1;
      	potq[jxyz].im *= ist.ngrid_1;

        fprintf(pf, "%.12f % .12f % .12f\n", sqrt(x2+y2+z2), potq[jxyz].re, potq[jxyz].im);
      }
    }
  }
  fclose(pf);

  free(potr);  free(kx2); free(ky2); free(kz2);

  return;
}

/************************************************************************/

#define SQRTPI       (sqrt(3.14159265358979323846))

double screenedcoulomb(double dr,double gamma)
{
  /*if (dr < EPSR) return (gamma);*/
  if (dr < EPSR) return (2.0*gamma/SQRTPI);
  return (erf(gamma * dr) / dr);
  /*return ((1.0 - exp(-gamma * dr)) / dr);*/
}

/************************************************************************/

void init_psi(zomplex *psi,lng_st ist,par_st par,long *idum)
{
  long jx, jy, jz, jzy, jxyz;
  long tidum = (*idum);

  for (jz = 0; jz < ist.nz; jz++) for (jy = 0; jy < ist.ny; jy++){
    for (jzy = ist.nx * (ist.ny * jz + jy), jx = 0; jx < ist.nx; jx++){
      jxyz = jzy + jx;
      psi[jxyz].re = (-1.0 + 2.0 * ran_nrc(&tidum));
      psi[jxyz].im = 0.0;
    }
  }
  normalize(psi,par.dv,ist.ngrid);
  (*idum) = tidum;
  return;
}

/************************************************************************/

double expot(double r, pot_st ppar) {
	return (ppar.e*(1.0/(exp((r-ppar.x0)/ppar.t)+1.0)-1.0));
}

/************************************************************************/
