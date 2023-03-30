/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "ar.h"
#include "unistd.h"

/*****************************************************************************/

void init_size(par_st *par, lng_st *ist) {
  FILE *pf; 
  long i, j, ieof, lenEvalParFile, count;
  double intExcEnergy, qpEnergy, sigma, *eval, *sige, a, holeBandEdgeE;
  char field[100], tmp[100];
  
  // Set defaults
  ist->nPseudoPot = 8192;          // size of pseudopotential files 
  ist->nStatesPerFilter = 8;       // number of states per filter    
  ist->nAtomTypes = 15;            // maximum number of atom types
  ist->nStochHotHoles = 40;        // sample 50 hot holes if stochastic = 1
  ist->nStochHotElecs = 40;        // sample 50 hot elecs if stochastic = 1
  ist->nStochOrbitals = 1600;      // number of stochastic orbitals used to approximate the Coulomb operator
  ist->nConsWindows = 10;          // default number of energy conservation windows
  ist->nMaxIntExcitons = 100;      // maximum number of interacting excitons that are used to make the biexcitonic states
  ist->seed = -234309823;          // random seed -> use this for testing
  ist->readInPsiEvalFiltFiles = 0; // default to begin completely new calculation. = 1 for already calculated psi-filt eval-filt files
  ist->readInHotQPs = 0;           // default to begin completely new calculation. = 1 for already calculated psiai and evalai
  ist->readAaiMatrices = 0;        // default to begin completely new calculation. = 1 for partly calculated AaiHole and AaiElec
  ist->nTrappedStates = 0;
  ist->debug = 0;                  // default to not printing/ running extra tests. debug = 1 -> verbose
  par->Elmin = -1;                 // minimum energy 
  par->Elmax = 1;                  // maximum energy
  par->Ekinmax = 10;               // kinetic energy cutoff 
  par->maxDeltaE = 0.0037;         // ~100 meV
  par->sigmaCutoff = 0.01;         // sigma=sqrt(<E^2>-<E>^2) -> values lower than sigmaCutoff are eigenstates of H
  par->temp = 298.0;               // electronic temperature
  par->fermiEnergy = -0.175;       // fermi energy, eigenstates with lower (higher) energies are holes (electrons) 

  // Read input.par if it exists - exit program otherwise 
  if ( access("input.par", F_OK) != -1 ) {
    pf = fopen("input.par", "r");
    i = 0;
    while (fscanf(pf, "%s", field) != EOF && i < 19) {
      if (! strcmp(field, "calcType")) fscanf(pf, "%s %s", tmp, &(ist->calcType));
      else if (! strcmp(field, "gridSize")) fscanf(pf, "%s %ld %ld %ld", tmp, &(ist->nx), &(ist->ny), &(ist->nz));
      else if (! strcmp(field, "nFilters"))  fscanf(pf, "%s %ld", tmp, &(ist->nFilterCycles)); 
      else if (! strcmp(field, "nNewtonIntSteps"))  fscanf(pf, "%s %ld", tmp, &(ist->nNewtonIntSteps)); 
      else if (! strcmp(field, "nStatesPerFilter"))  fscanf(pf, "%s %ld", tmp, &(ist->nStatesPerFilter));
      else if (! strcmp(field, "nThreads"))  fscanf(pf, "%s %ld", tmp, &(ist->nThreads)); 
      else if (! strcmp(field, "seed"))  fscanf(pf, "%s %ld", tmp, &(ist->seed));
      else if (! strcmp(field, "nTrappedStates"))  fscanf(pf, "%s %ld", tmp, &(ist->nTrappedStates)); 
      else if (! strcmp(field, "nStochOrbitals"))  fscanf(pf, "%s %ld", tmp, &(ist->nStochOrbitals));
      else if (! strcmp(field, "nStochHotHoles"))  fscanf(pf, "%s %ld", tmp, &(ist->nStochHotHoles));
      else if (! strcmp(field, "nStochHotElecs"))  fscanf(pf, "%s %ld", tmp, &(ist->nStochHotElecs));      
      else if (! strcmp(field, "nConsWindows"))  fscanf(pf, "%s %ld", tmp, &(ist->nConsWindows));
      else if (! strcmp(field, "maxDeltaE"))  fscanf(pf, "%s %lg", tmp, &(par->maxDeltaE)); 
      else if (! strcmp(field, "sigmaCutoff"))  fscanf(pf, "%s %lg", tmp, &(par->sigmaCutoff));
      else if (! strcmp(field, "temp"))  fscanf(pf, "%s %lg", tmp, &(par->temp));
      else if (! strcmp(field, "fermiEnergy"))  fscanf(pf, "%s %lg", tmp, &(par->fermiEnergy));
      else if (! strcmp(field, "nMaxIntExcitons")) fscanf(pf, "%s %ld", tmp, &(ist->nMaxIntExcitons));
      else if (! strcmp(field, "readInPsiEvalFiltFiles")) fscanf(pf, "%s %ld", tmp, &(ist->readInPsiEvalFiltFiles));
      else if (! strcmp(field, "readInHotQPs")) fscanf(pf, "%s %ld", tmp, &(ist->readInHotQPs));
      else if (! strcmp(field, "readAaiMatrices")) fscanf(pf, "%s %ld", tmp, &(ist->readAaiMatrices)); 
      else if (! strcmp(field, "debug"))  fscanf(pf, "%s %ld", tmp, &(ist->debug));
      else {
        printf("Invalid input field and/or format - equal sign required after each field\n");
        printf("Only allowed fields are (case-sensitive):\n\n");
        printf("calcType = sNI (required, dNI, sNI, testNI, dI, s1I, test1I)\n");
        printf("gridSize = nx ny nz (required, all integers)\n");
        printf("nFilters = 32 (required, integer, should be multiple of nThreads)\n");
        printf("nNewtonIntSteps = 2048 (required, integer, normally multiple of 2)\n");
        printf("nStatesPerFilter = 8 (optional, 8 default, normally 8-20)\n");
        printf("nThreads = 1 (optional, serial default, number of openmp threads)\n");
        printf("seed = -234309823 (optional, integer, random number)\n");
        printf("nStochOrbitals = 1600 (optional, 1600 default, number of stochastic orbitals in ~Coulomb operator)\n");
        printf("nStochHotHoles = 50 (optional, 50 default, number of stochastic hot hole orbitals)\n");
        printf("nStochHotElecs = 50 (optional, 50 default, number of stochastic hot elec orbitals)\n");
        printf("nConsWindows = 10 (optional, 10 default, number of energy conservation windows)\n");
        printf("maxDeltaE = 0.0037 (optional, ~100 meV default, largest energy conservation window)\n");
        printf("sigmaCutoff = 0.01 (optional, 0.01 default, criteria to determine if eigenstate or not)\n");
        printf("temp = 298 (optional, 298 K default, electronic temperature in final AR lifetime calculation)\n");
        printf("fermiEnergy = -0.175 (optional, -0.175 default, eigenstates with lower (higher) energies are holes (electrons)\n");
        printf("nMaxIntExcitons = 100 (optional, 100 default, max number of excitons used to make the initial biexcitonic states\n");
        printf("readInPsiEvalFiltFiles = 1 (optional, 0 is default, reads in psi-filt- to pick up where filtering left off\n");
        printf("readInHotQPs = 1 (optional, 0 is default, reads in psiaiEigOnly and evalaiEigOnly files when = 1)\n");
        printf("readAaiMatrices = 1 (optional, 0 is default, reads in AaiMatrices.dat file when = 1)\n");
        printf("debug = 1 (optional, 0 is default, prints out more information when debug = 1)\n");
        fflush(stdout);
        exit(EXIT_FAILURE);
      }
      i++;
    }
    fclose(pf);
  }
  else {
    printf("\n\nNo input.par file detected in current working directory - the program is exiting!!!\n\n");
    exit(EXIT_FAILURE);
  }

  // Set true false parameters based on calcType
  setCalcTypeParameters(ist);

  // Read the number of atoms from conf.par if it exits - exit program otherwise
  if ( access("conf.par", F_OK) != -1 ) {
    pf = fopen("conf.par", "r");
    fscanf(pf, "%ld", &ist->nAtoms); 
    fclose(pf);
  }
  else {
    printf("\n\nNo conf.par file detected in current working directory - the program is exiting!!!\n\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
  }

  // Read in eval.par and determine number of eigenstates if it exits - exit program otherwise
  ist->totalHomo = ist->totalLumo = 0;
  if ( access("eval.par", F_OK) != -1 ) {
    pf = fopen("eval.par" , "r");
    for (i = ieof = 0; ieof != EOF; i++) {
      ieof = fscanf(pf, "%ld %lg %lg", &j, &qpEnergy, &sigma);
      if (sigma < par->sigmaCutoff && qpEnergy < par->fermiEnergy) {
        par->homoEnergy = qpEnergy;
        ist->homoIndex = i;
        ist->totalHomo++;
      }
      else if (sigma < par->sigmaCutoff && qpEnergy > par->fermiEnergy && ieof != EOF) {
        if (! ist->totalLumo) {
          par->lumoEnergy = qpEnergy;
          ist->lumoIndex = i;  
        }
        ist->totalLumo++;
      } 
    }
    fclose(pf);
    lenEvalParFile = i-1; // length of eval.par
  }
  else {
    printf("\n\nNo eval.par file detected in current working directory - the program is exiting!!!\n\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
  }
  
  // Determine number of hole and electron states for which Coulomb matrix elements will need to be calculated
  par->kbT = KB*par->temp/AUTOEV;
  par->boltzEnergyRange = 3.0*par->kbT;
  ist->nHoles = ist->nElecs = 0;
  if (ist->intAR) {
    ist->nHoles = ist->totalHomo;
    ist->nElecs = ist->totalLumo;
  }
  else if (ist->nonIntAR) {
    eval = (double *) calloc(lenEvalParFile, sizeof(double)); 
    sige = (double *) calloc(lenEvalParFile, sizeof(double)); 
    pf = fopen("eval.par" , "r");
    for (i = 0; i < lenEvalParFile; i++){
      fscanf(pf, "%ld %lg %lg", &j, &eval[i], &sige[i]);
    }
    fclose(pf);
    // count hole states
    count = 0;
    holeBandEdgeE = par->homoEnergy;
    for (i = ist->homoIndex; i >= 0; i--) {
      if ((eval[i] >= (holeBandEdgeE-par->boltzEnergyRange)) && sige[i] < par->sigmaCutoff) {
        ist->nHoles++;
      }
      // get band edge hole state
      if (sige[i] < par->sigmaCutoff) {
          if (count == ist->nTrappedStates) {
              holeBandEdgeE = eval[i];
              if (ist->nTrappedStates > 0) {
                ist->nHoles++;
              }
          }
          count++;
      }
    }
    // count electron states
    for (i = ist->lumoIndex; i < lenEvalParFile; i++) {
      if ((eval[i] <= par->lumoEnergy+par->boltzEnergyRange) && sige[i] < par->sigmaCutoff) {
        ist->nElecs++;
      }
    }

    free(eval);  free(sige);
  }

  // Determine the number of hot hole and electron states if readInHotQPs = 1 (true)
  if (ist->readInHotQPs) {
    ist->nHotHoles = 0; ist->nHotElecs = 0;
    if ( access("evalai.par", F_OK) != -1 ) {
      pf = fopen("evalai.par" , "r");
      for (i = ieof = 0; ieof != EOF; i++) {
        ieof = fscanf(pf, "%ld %lg %lg", &j, &qpEnergy, &sigma);
        if (sigma < par->sigmaCutoff && qpEnergy < par->fermiEnergy) {
          ist->nHotHoles++;
        }
        else if (sigma < par->sigmaCutoff && qpEnergy > par->fermiEnergy && ieof != EOF) {
          ist->nHotElecs++;
        } 
      }
      fclose(pf);
      ist->itot = ist->nHotHoles;
      ist->atot = ist->nHotElecs;
      ist->nHotNonIntStates = ist->nHotHoles + ist->nHotElecs;
    }
    else {
      printf("\n\nNo evalai.par file detected in current working directory - the program is exiting!!!\n\n");
      fflush(stdout);
      exit(EXIT_FAILURE);
    }    
  }

  // Set useful parameters based on input parameters
  ist->nHolesPlusElecs = ist->nHoles+ist->nElecs;
  ist->nNonIntExcitons = ist->nHoles*ist->nElecs; 
  ist->nFilteredStates = ist->nStatesPerFilter * ist->nFilterCycles;  // total number of filtered states
  ist->nx_1 = 1.0 / (double)(ist->nx);
  ist->ny_1 = 1.0 / (double)(ist->ny);
  ist->nz_1 = 1.0 / (double)(ist->nz);
  ist->nGridPoints = ist->nx * ist->ny * ist->nz; 
  ist->nGridPoints_1 = 1.0 / (double)(ist->nGridPoints);
  par->fundamentalGap = par->lumoEnergy - par->homoEnergy;
  if (ist->intAR) {
    ist->nIntExcitons = 0;
    if ( access("exciton.par", F_OK) != -1 ) {
      pf = fopen("exciton.par", "r");
      for (i = ieof = 0; ieof != EOF; i++) {
  	ieof = fscanf(pf, "%ld %lg %lg %lg %lg", &j, &intExcEnergy, &a, &a, &a);
        if (! i) {
    	  par->opticalGap = intExcEnergy;
    	  par->minInitE = 2.0*par->opticalGap;
    	  par->maxInitE = par->minInitE + par->boltzEnergyRange;
    	  ist->nIntExcitons++;
          // printf("%i %.8lf %i l0\n", i, intExcEnergy, ist->nIntExcitons);
          // fflush(stdout);
    	}
        else if (i < ist->nTrappedStates) {
          // DJ edits:
          ist->nIntExcitons++;
          // if (i == 1) {
          //   par->minInitE = par->opticalGap + intExcEnergy;
          // }
          // par->maxInitE = par->opticalGap + intExcEnergy + par->boltzEnergyRange;
          // ist->nIntExcitons++;
          // printf("%i %.8lf %i l1\n", i, intExcEnergy, ist->nIntExcitons);
          // fflush(stdout);
        }
        // DJ
        else if (i == ist->nTrappedStates) {
          par->minInitE = par->opticalGap + intExcEnergy;
          par->maxInitE = par->minInitE + par->boltzEnergyRange;
          ist->nIntExcitons++;
          // printf("%i %.8lf %i %.8lf %.8lf l2\n", i, intExcEnergy, ist->nIntExcitons, par->minInitE, par->maxInitE);
          // fflush(stdout);
        }
        // DJ
    	// else if (intExcEnergy <= (par->opticalGap + par->boltzEnergyRange) && (ist->nIntExcitons < ist->nMaxIntExcitons)) {
    	else if (intExcEnergy <= (par->maxInitE - par->opticalGap) && (ist->nIntExcitons < ist->nMaxIntExcitons)) {
          ist->nIntExcitons++;
          // printf("%i %.8lf %i l3\n", i, intExcEnergy, ist->nIntExcitons);
          // fflush(stdout);
    	}
    	else {
          // printf("%i %.8lf %i breaking...\n", i, intExcEnergy, ist->nIntExcitons);
          // fflush(stdout);
    	  break; // all other excitonic states are of too high energy
    	}
      } 
      fclose(pf);
    }
    else {
      printf("\n\nNo exciton.par file detected in current working directory - the program is exiting!!!\n\n");
      exit(EXIT_FAILURE);
    }
  }
  else if (ist->nonIntAR) {
    // DJ edits
    if (ist->nTrappedStates > 0) {
        par->minInitE = par->fundamentalGap + par->lumoEnergy - holeBandEdgeE;
    }
    else{
        par->minInitE = 2.0*par->fundamentalGap;
    }
    par->maxInitE = par->minInitE + par->boltzEnergyRange;  
  }

  // Write the input parameters to stdout
  writeInputParameters(*ist, *par, stdout);
  fflush(stdout);

  return;
}

/*****************************************************************************/

void setCalcTypeParameters(lng_st *ist) {

  if (! strcmp(ist->calcType, "dNI")) {
    ist->deterministic = 1; ist->stochastic = 0; 
    ist->nonIntAR = 1; ist->intAR = 0;  
  }
  else if (! strcmp(ist->calcType, "sNI")) {
    ist->deterministic = 0; ist->stochastic = 1; 
    ist->nonIntAR = 1; ist->intAR = 0;  
  }
  else if (! strcmp(ist->calcType, "dI")) {
    ist->deterministic = 1; ist->stochastic = 0; 
    ist->nonIntAR = 0; ist->intAR = 1;
  }
  else if (! strcmp(ist->calcType, "s1I")) {
    ist->deterministic = 0; ist->stochastic = 1; 
    ist->nonIntAR = 0; ist->intAR = 1;
  }
  else if (! strcmp(ist->calcType, "s2I")) {
    ist->deterministic = 0; ist->stochastic = 2; 
    ist->nonIntAR = 0; ist->intAR = 1;
  }
  else if (! strcmp(ist->calcType, "s3I")) {
    ist->deterministic = 0; ist->stochastic = 3; 
    ist->nonIntAR = 0; ist->intAR = 1;
  }
  else if (! strcmp(ist->calcType, "testNI")) {
    ist->deterministic = 1; ist->stochastic = 1; 
    ist->nonIntAR = 1; ist->intAR = 0;
  }
  else if (! strcmp(ist->calcType, "test1I")) {
    ist->deterministic = 1; ist->stochastic = 1; 
    ist->nonIntAR = 0; ist->intAR = 1;
  }
  else if (! strcmp(ist->calcType, "test2I")) {
    ist->deterministic = 1; ist->stochastic = 2; 
    ist->nonIntAR = 0; ist->intAR = 1;
  }
  else if (! strcmp(ist->calcType, "test3I")) {
    ist->deterministic = 1; ist->stochastic = 3; 
    ist->nonIntAR = 0; ist->intAR = 1;
  }
  else {
    printf("An incorrect calcType was entered - the program is exiting!\n"); 
    exit(EXIT_FAILURE);
  }

  return;
}

/*****************************************************************************/

void init(double *vx, double *vy, double *vz, double *ksqr, double *potl, 
          double *rx, double *ry, double *rz, par_st *par, lng_st *ist)
{
  FILE *pf; atm_st *atm;
  long jexc, jx, jy, jz, jyz, jxyz, iatom, ntot, tmp, *nPseudoPot, flags=0;
  double del, mx, my, mz, xd, yd, zd, dx, dy, dz, *ksqrx, *ksqry, *ksqrz;
  double *vr, *potatom, *dr, sum, rex, potEx; 
  pot_st ppar;

  if ((ksqrx  = (double*)calloc(ist->nx,sizeof(double)))==NULL)nerror("ksqrx");
  if ((ksqry  = (double*)calloc(ist->ny,sizeof(double)))==NULL)nerror("ksqry");
  if ((ksqrz  = (double*)calloc(ist->nz,sizeof(double)))==NULL)nerror("ksqrz");
  if ((atm = (atm_st*)calloc(ist->nAtoms,sizeof(atm_st)))==NULL)nerror("atm");
  
  // Read in the external potential and store it in ppar if it exists
  if ( access( "expot.par", F_OK) != -1 ) {
    pf = fopen("expot.par", "r");
    fscanf(pf, "%lg", &ppar.e);  /*** external potential parameter: Energy Scale ***/
    fscanf(pf, "%lg", &ppar.x0); /*** external potential parameter: Ef ***/
    fscanf(pf, "%lg", &ppar.t);  /*** external potential parameter: kT ***/
    fclose(pf);
    printf("expot.par was detected\n");
    printf("External potential energy scale, e = % .6f\n", ppar.e);
    printf("External potential Fermi radius, x0 = % .6f\n", ppar.x0);
    printf("External potential steepness (Temp), kT = % .6f\n", ppar.t);
  } else {
    ppar.e = ppar.x0 = 0.0; // setting ppar.e = 0 results in ext pot being 0
    ppar.t = 1.0;
    printf("No expot.par file detected in cwd - setting external potential to 0\n");
  }

  // Read the passivated nanocrystal configuration 
  pf = fopen("conf.par" , "r");
  fscanf(pf, "%ld", &tmp);
  ist->mfermi = read_conf(rx, ry, rz, atm, ist->nAtoms, pf);
  fclose(pf);

  // Determine the box size
  xd = rint(0.5 * get_dot_ligand_size_z(rx, ist->nAtoms) + 5.0);
  yd = rint(0.5 * get_dot_ligand_size_z(ry, ist->nAtoms) + 5.0);
  zd = rint(0.5 * get_dot_ligand_size_z(rz, ist->nAtoms) + 5.0);
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

  /***initializing the potential vector  ***/
  for (jx = 0, dx = par->xmin; jx < ist->nx; jx++, dx += par->dx) vx[jx] = dx;
  for (jy = 0, dy = par->ymin; jy < ist->ny; jy++, dy += par->dy) vy[jy] = dy;
  for (jz = 0, dz = par->zmin; jz < ist->nz; jz++, dz += par->dz) vz[jz] = dz;
  
  /*** read pseudopotentials ***/
  dr  = (double *) calloc(ist->nAtomTypes, sizeof(double));
  vr  = (double *) calloc(ist->nPseudoPot*ist->nAtomTypes, sizeof(double));
  potatom = (double *) calloc(ist->nPseudoPot*ist->nAtomTypes, sizeof(double));
  nPseudoPot = (long *) calloc(ist->nAtomTypes, sizeof(long));
  read_pot(vr, potatom, nPseudoPot, dr, atm, ist->nPseudoPot, ist->nAtomTypes);

  par->Vmin = 1.0e10;
  par->Vmax = -1.0e10;

  omp_set_dynamic(0);
  omp_set_num_threads(ist->nThreads);
#pragma omp parallel for private(dx,dy,dz,del,jy,jx,jyz,jxyz,sum,iatom)
  for (jz = 0; jz < ist->nz; jz++) {
    for (jy = 0; jy < ist->ny; jy++) {
      jyz = ist->nx * (ist->ny * jz + jy);
      for (jx = 0; jx < ist->nx; jx++) {
        rex = sqrt(vx[jx]*vx[jx]+vy[jy]*vy[jy]+vz[jz]*vz[jz]);
        potEx = expot(rex, ppar);
      	jxyz = jyz + jx;
      	for (sum = 0.0, iatom = 0; iatom < ist->nAtoms; iatom++) {
      	  dx = vx[jx] - rx[iatom];
      	  dy = vy[jy] - ry[iatom];
      	  dz = vz[jz] - rz[iatom];
      	  del = sqrt(dx * dx + dy * dy + dz * dz);
      	  sum += interpolate(del,dr[atm[iatom].natyp],vr,potatom,ist->nPseudoPot,nPseudoPot[atm[iatom].natyp],atm[iatom].natyp);
      	}
      	potl[jxyz] = sum+potEx;
      }
    }
  }

  // Print out the external potential as a function of the absolute value of the z coordinate
  pf = fopen("pS.dat", "w");
  for (jz = 0; jz < ist->nz; jz++) {
    potEx = expot(fabs(vz[jz]), ppar);
    fprintf(pf, "% .6f % .8f\n", vz[jz], potEx);
  }
  fclose(pf);  

  // Get max and min values of the potential
  for (jxyz= 0; jxyz < ist->nGridPoints; jxyz++) {
    if (par->Vmax < potl[jxyz]) par->Vmax = potl[jxyz];
    if (par->Vmin > potl[jxyz]) par->Vmin = potl[jxyz];
  }
  printf("Potential energy statistics: dV = %.5f Vmin = % .5f Vmax = %.5f\n",
	  par->Vmax-par->Vmin, par->Vmin, par->Vmax);

  par->dE = (0.5 * sqr(PIE) / (mx*par->dx*par->dx) +
             0.5 * sqr(PIE) / (my*par->dy*par->dy) +    
             0.5 * sqr(PIE) / (mz*par->dz*par->dz));
  printf("dT = %g\n", par->dE);
  fflush(stdout);

  // Free dynamically allocated memory
  free(dr); free(vr); free (atm);  free(nPseudoPot);

  return;
}

/****************************************************************************/

void init_pot(double *vx, double *vy, double *vz, zomplex *potq, par_st par, lng_st ist,
              fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi)
{
  FILE *pf;
  long jx, jy, jz, jyz, jxyz, sx, sy, sz;
  double dr, x2, y2, z2, *kx2, *ky2, *kz2, alpha, cosa, sina;
  double boxl = (double)(ist.nx) * par.dx;
  zomplex *potr, tmp;

  // Dynamically allocate memory
  if ((kx2  = (double *) calloc(ist.nx, sizeof(double))) == NULL) nerror("kx2");
  if ((ky2  = (double *) calloc(ist.ny, sizeof(double))) == NULL) nerror("ky2");
  if ((kz2  = (double *) calloc(ist.nz, sizeof(double))) == NULL) nerror("kz2");
  if ((potr = (zomplex *) calloc(ist.nGridPoints, sizeof(zomplex))) == NULL) nerror("potr");
  
  // Calculate the real space Coulomb potential
  for (jx = 0; jx < ist.nGridPoints; jx++) {
    potr[jx].re = potr[jx].im = 0.0;
  }
  for (jz = 0; jz < ist.nz; jz++) {
    z2 = sqr(vz[jz]);
    for (jy = 0; jy < ist.ny; jy++) {
      y2 = sqr(vy[jy]);
      jyz = ist.nx * (ist.ny * jz + jy);
      for (jx = 0; jx < ist.nx; jx++) {
      	x2 = sqr(vx[jx]);
      	jxyz = jyz + jx;
      	dr = sqrt(x2 + y2 + z2);
      	if (dr < boxl) potr[jxyz].re = screenedcoulomb(dr, par.gamma);
      }
    }
  }

  // Fourier transform the r-space potential 
  memcpy(&fftwpsi[0], &potr[0], ist.nGridPoints*sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  memcpy(&potq[0], &fftwpsi[0], ist.nGridPoints*sizeof(potq[0]));
  
  for (kx2[0] = 0.0, jx = 1; jx <= ist.nx / 2; jx++) {
    kx2[jx] = (kx2[ist.nx-jx] = sqr((double)(jx) * par.dkx));
  }
  for (ky2[0] = 0.0, jy = 1; jy <= ist.ny / 2; jy++) {
    ky2[jy] = (ky2[ist.ny-jy] = sqr((double)(jy) * par.dky));
  }
  for (kz2[0] = 0.0, jz = 1; jz <= ist.nz / 2; jz++) {
    kz2[jz] = (kz2[ist.nz-jz] = sqr((double)(jz) * par.dkz));
  }

  //pf = fopen("qSpaceCoulombPot.dat", "w");
  for (jz = 0; jz < ist.nz; jz++) {
    z2 = kz2[jz];
    for (jy = 0; jy < ist.ny; jy++) {
      y2 = ky2[jy];
      jyz = ist.nx * (ist.ny * jz + jy);
      for (jx = 0; jx < ist.nx; jx++) {
      	x2 = kx2[jx];
      	jxyz = jyz + jx;
      	alpha = PIE * (double)(jx + jy + jz + ist.nGridPoints / 2);
      	cosa = cos(alpha);
      	sina = sin(alpha);
	
       	tmp.re = potq[jxyz].re;
      	tmp.im = potq[jxyz].im;

      	potq[jxyz].re = (tmp.re * cosa - tmp.im * sina) * par.dv;
      	potq[jxyz].im = (tmp.re * sina + tmp.im * cosa) * par.dv;
      	
        // See equation (4.3) in Martyna & Tuckerman. JCP, 110, 2810 (1999)
        // Full paper: 99_JCP_Tuckerman_KSpaceMethodLongRangeInteractions.pdf
        // W(g) = [4*pi/g^2]*exp(-g^2/(4*alphaewd^2)) + phiScreenedCoulomb(g)
        double  q2 = x2 + y2 + z2 + 1.0e-20;
      	potq[jxyz].re += (FOURPI * (1.0 - exp(-0.25* q2 / par.gamma2)) / q2);


      	potq[jxyz].re *= ist.nGridPoints_1;
      	potq[jxyz].im *= ist.nGridPoints_1;

      	//fprintf(pf, "%.12f % .12f % .12f\n", sqrt(x2+y2+z2), potq[jxyz].re, potq[jxyz].im);
      }
    }
  }
  //fclose(pf);

  // Debug printing
  // writeComplexArray(potr, ist.nGridPoints, "rSpaceCoulombPot");

  // Free dynamically allocated memory
  free(potr);  free(kx2); free(ky2); free(kz2);

  return;
}

/************************************************************************/

#define SQRTPI       (sqrt(3.14159265358979323846))

double screenedcoulomb(double dr, double gamma)
{
  /*if (dr < EPSR) return (gamma);*/
  if (dr < EPSR) return (2.0*gamma/SQRTPI);
  return (erf(gamma * dr) / dr);
  /*return ((1.0 - exp(-gamma * dr)) / dr);*/
}

/************************************************************************/

void init_psi(zomplex *psi, lng_st ist, par_st par, long *idum)
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
  normalize(psi,par.dv,ist.nGridPoints);
  (*idum) = tidum;
  return;
}

/************************************************************************/

double expot(double r, pot_st ppar) {
  return (ppar.e*(1.0/(exp((r-ppar.x0)/ppar.t)+1.0)-1.0));
}

/************************************************************************/
