/*****************************************************************************/

#include "fd.h"
#include <float.h>

/*****************************************************************************/

int main(int argc, char *argv[]) {
    FILE *ppsi, *peval;  
    long i, a, j, thomo, tlumo, indexfirsthomo, flags=0;
    double *eval, *de, *ksqr, *vx, *vy, *vz, *psi, *psidummy, *rx, *ry, *rz, *poth;
    double *potl, *bsmat, *h0mat, *mux, *muy, *muz, *mx, *my, *mz, *rs;
    zomplex *potq, *potqx;
    fftw_plan_loc *planfw, *planbw; fftw_complex *fftwpsi;
    par_st par;
    long_st ist;

    /*************************************************************************/
    writeCurrentTime(stdout);
    writeSeparation(stdout);
    init_size(argc, argv, &par, &ist);
  
    /*************************************************************************/
    fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.ngrid*ist.nthreads);
    potl  = (double *) calloc(ist.ngrid, sizeof(double));
    potq  = (zomplex *) calloc(ist.ngrid, sizeof(zomplex));
    potqx  = (zomplex *) calloc(ist.ngrid, sizeof(zomplex));
    poth = (double *) calloc(ist.ngrid*ist.nthreads, sizeof(double));
    ksqr = (double *) calloc(ist.ngrid, sizeof(double));
    vx = (double *) calloc(ist.nx, sizeof(double));
    vy = (double *) calloc(ist.ny, sizeof(double));
    vz = (double *) calloc(ist.nz, sizeof(double));
    rx = (double *) calloc(ist.natom, sizeof(double));
    ry = (double *) calloc(ist.natom, sizeof(double));
    rz = (double *) calloc(ist.natom, sizeof(double));
  
    /**************************************************************************/
    init(potl, vx, vy, vz, ksqr, rx, ry, rz, &par, &ist);

    /*** initialization for the fast Fourier transform ***/
    fftw_plan_with_nthreads(ist.nthreads);
  
    planfw = (fftw_plan_loc *) calloc(ist.nthreads, sizeof(fftw_plan_loc));
    planbw = (fftw_plan_loc *) calloc(ist.nthreads, sizeof(fftw_plan_loc));
    for (i = 0; i < ist.nthreads; i++) { 
        planfw[i] = fftw_plan_dft_3d(ist.nz, ist.ny, ist.nx, &fftwpsi[i*ist.ngrid], 
                                     &fftwpsi[i*ist.ngrid], FFTW_FORWARD, flags);
        planbw[i] = fftw_plan_dft_3d(ist.nz, ist.ny, ist.nx, &fftwpsi[i*ist.ngrid],
                                     &fftwpsi[i*ist.ngrid], FFTW_BACKWARD, flags);
    }
    init_pot(vx, vy, vz, potq, potqx, par, ist, planfw[0], planbw[0], &fftwpsi[0]);

    /*************************************************************************/
    eval = (double *) calloc(ist.nhomo+1, sizeof(double)); 
    de = (double *) calloc(ist.nhomo+1, sizeof(double)); 

    peval = fopen("eval.par" , "r");
    for (i = 0; i < ist.nhomo+1; i++)
        fscanf (peval,"%ld %lg %lg",&a,&eval[i],&de[i]);
    fclose (peval);

    peval = fopen("eval.dat" , "w");
    for (i = 0; i < ist.nhomo+1; i++){
        fprintf (peval,"%ld %g %g\n",i,eval[i],de[i]);
    }
    fclose(peval);
  
    for (thomo = 0, i = ist.nhomo; i >= 0; i--){
        if (de[i] < par.deps) thomo++;
        if (thomo == ist.totalhomo) {
            indexfirsthomo = i;
            break;
        }
    }

    printf("The index of lowest energy occupied level used = %ld\n", indexfirsthomo); 

    free(eval);
    free(de);

    psidummy = (double*)calloc(ist.ngrid,sizeof(double));
    eval = (double*)calloc(ist.ms,sizeof(double)); 
    de = (double*)calloc(ist.ms,sizeof(double)); 

    /**********************************************************************/
    /*              CHANGES TO THE FOLLOWING CODE START HERE              */ 
    /**********************************************************************/
    /* Ideally I would prefer it to be outside of main, as its logic is more 
     * specific than "run the BSE Calc", but in the interest of time I'm 
     * going to leave it like this and focus back on calculating auger rates
     */

    /**********************************************************************/
    /*               Read in eval.par and store it in memory              */
    /**********************************************************************/
    struct eindex { // Struct I used to store index evals and deps in memory
        long index;
        double evalue;
        double deps;
    }; 
    
    struct eindex *evalindex;
    int size = 1024; // Initial size of array in bytes, is allowed to resize
    long counter = 0;

    if ((evalindex = calloc(size, sizeof(struct eindex))) == NULL) {
        fprintf(stderr, "Memory error in main\n");
        exit(EXIT_FAILURE);
    }

    peval = fopen("eval.par" , "r");
    if (peval) {
         while (fscanf(peval, "%ld %lg %lg", &evalindex[counter].index, 
                &evalindex[counter].evalue, &evalindex[counter].deps) == 3) {
            counter++;
            if (counter == size - 1) {
                size *= 2;
                // TODO: This does not check if realloc succeeded or not
                evalindex = realloc(evalindex, size * sizeof(struct eindex));
            }
        }
    }
    fclose(peval);

    // Read ead in psi.par as we have already read in psi.par
    ppsi = fopen("psi.par" , "r");

    double *psihomo = calloc(ist.nhomo*ist.ngrid, sizeof(double));
    double *psilumo = calloc(ist.nlumo*ist.ngrid, sizeof(double));
    if (!psihomo || !psilumo) terminate("Failed to allocate memory for psihomo/psilumo");

    long foffset = ist.ngrid * sizeof(double);  // for random access 
    char fname[80] = {0};
    long nstates = 0;  // Total number of states

    fseek(ppsi, foffset * ist.nhomo, SEEK_SET); 
    counter = ist.nhomo; // Set loop counter to HOMO index
    thomo = 0;           // Counter for hole states used
    while (counter >= indexfirsthomo && thomo < ist.maxHoleStates) {
        fread(psidummy, sizeof(double), ist.ngrid, ppsi);
        if (evalindex[counter].deps < par.deps) {    
            normalize(psidummy, par.dv, ist.ngrid);
            sprintf(fname, "pzv%d.dat", counter);
			if ((z_project(psidummy, vz, par, ist, fname) == 0)) {
                eval[nstates] = evalindex[counter].evalue;
                de[nstates] = evalindex[counter].deps;
                for (j = 0; j < ist.ngrid; j++) {
                    psihomo[thomo * ist.ngrid + j] = psidummy[j];
                }
                nstates++;
                thomo++;
			}
        }
        counter--;
        fseek(ppsi, -2 * foffset, SEEK_CUR); 
    }

    fseek(ppsi, foffset * (ist.nhomo + 1), SEEK_SET);
    counter = ist.nhomo + 1;
    tlumo = 0;

    while (tlumo < ist.maxElecStates && tlumo != ist.totallumo) {
        fread(psidummy, sizeof(double), ist.ngrid, ppsi);
        if (evalindex[counter].deps < par.deps) {
            normalize(psidummy, par.dv, ist.ngrid);
            sprintf(fname, "pzc%ld.dat", counter);
            if ((z_project(psidummy, vz, par, ist, fname)) == 0) {
			    eval[nstates] = evalindex[counter].evalue;
                de[nstates] = evalindex[counter].deps;
                for (j = 0; j < ist.ngrid; j++) {      
                    psilumo[tlumo * ist.ngrid + j] = psidummy[j];
                }
            	nstates++;
            	tlumo++;
			}
        }
        counter++;
    }
    psi = (double *) calloc(nstates, ist.ngrid*sizeof(double));
    if (!psi) terminate("Failed to allocate memory for psi");

    /**********************************************************************/

    for (i = thomo - 1, a = 0; i >= 0, a < thomo; i--, a++) {
        for (j = 0; j < ist.ngrid; j++) {
            psi[a * ist.ngrid + j] = psihomo[i * ist.ngrid + j];
        }
    }

    for (i = thomo, a = 0; i < tlumo + thomo, a < tlumo; i++, a++) {
        for (j = 0; j < ist.ngrid; j++) {
            psi[i * ist.ngrid + j] = psilumo[a * ist.ngrid + j];
        }
    }

    double tmp1, tmp2;
    /* Rearrange eval */
    for (i = 0, j = thomo - 1; i < j; i++, j--) {
        tmp1 = eval[i];    tmp2 = de[i];
        eval[i] = eval[j]; de[i] = de[j];
        eval[j] = tmp1;    de[j] = tmp2;
    }
    /* Write smaller psi.par and eval.par for augerBSE */
    FILE *new_eval = fopen("BSEeval.par", "w");
    FILE *new_psi = fopen("BSEpsi.par", "w");
    if (!new_eval || !new_psi) terminate("Failed opening BSEeval.par");

    fwrite(psi, sizeof(double) * ist.ngrid, nstates, new_psi);

    for (i = 0; i < nstates; i++) {
        fprintf(new_eval, "%ld %.*g %.*g\n", i, DBL_DIG, eval[i], DBL_DIG, de[i]); 
    }
    fclose(new_eval);
    fclose(new_psi);

    /**********************************************************************/
    fclose(ppsi);
    free(psidummy);
    free(evalindex);
    free(psihomo);
    free(psilumo);
  
    printf("Number of hole eigenstates used in the BSE calculation = %ld\n", thomo);
    printf("Number of electron eigenstates used in the BSE calculation =  %ld\n", tlumo);
    printf("Total number of eigenstates used in the BSE calculation = %ld\n", nstates);

    /*************************************************************************/
    ist.ms = nstates;
    ist.nlumo = thomo;
    ist.nhomo = thomo - 1;
    ist.totallumo = tlumo;
    ist.totalhomo = thomo;
  
    /*                              END CHANGES                           */
    /*************************************************************************/

    /**************************************************************************/
    //get_energy_range(vx, vy, vz, ksqr, potl, &par, ist, planfw, planbw, fftwpsi);
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
         vjbai and vabji are the coulomb matrix elements need to be used to
         generate the spin-depedent matrix elements as described by
         the last equation in our codument.  ***/

    ist.ms2 = ist.totallumo * ist.totalhomo;
    bsmat = (double *) calloc(ist.ms2*ist.ms2, sizeof(double)); 
    h0mat = (double *) calloc(ist.ms2*ist.ms2, sizeof(double)); 
    
    single_coulomb_openmp(psi, potq, potqx, poth, eval, ist, par, planfw, planbw, fftwpsi, bsmat, h0mat);

    ppsi = fopen("bs.dat", "w");
    for (i = 0; i < ist.ms2; i++, fprintf(ppsi,"\n"))
        for (j = 0; j < ist.ms2; j++)
            fprintf(ppsi,"%.*g ", DBL_DIG, bsmat[i*ist.ms2+j]);
    fclose(ppsi);

    ppsi = fopen("h0.dat", "w");
    for (i = 0; i < ist.ms2; i++, fprintf(ppsi,"\n"))
        for (j = 0; j < ist.ms2; j++)
             fprintf(ppsi,"%.*g ", DBL_DIG, h0mat[i*ist.ms2+j]);
    fclose(ppsi);

    mux = (double *) calloc(ist.totallumo*ist.totalhomo, sizeof(double)); // <psi_i|ux|psi_a>
    muy = (double *) calloc(ist.totallumo*ist.totalhomo, sizeof(double)); // <psi_i|uy|psi_a>
    muz = (double *) calloc(ist.totallumo*ist.totalhomo, sizeof(double)); // <psi_i|uz|psi_a>
    mx  = (double *) calloc(ist.totallumo*ist.totalhomo, sizeof(double)); // <psi_i|mx|psi_a>
    my  = (double *) calloc(ist.totallumo*ist.totalhomo, sizeof(double)); // <psi_i|my|psi_a>
    mz  = (double *) calloc(ist.totallumo*ist.totalhomo, sizeof(double)); // <psi_i|mz|psi_a>
    rs  = (double *) calloc(ist.totallumo*ist.totalhomo, sizeof(double)); // <psi_a|u|psi_i>.<psi_a|m|psi_i>
    // px_ni  = (double *) calloc(ist.totallumo+ist.totalhomo, sizeof(double)); // <psi|px|psi>
    // py_ni  = (double *) calloc(ist.totallumo+ist.totalhomo, sizeof(double)); // <psi|py|psi>
    // pz_ni  = (double *) calloc(ist.totallumo+ist.totalhomo, sizeof(double)); // <psi|pz|psi>

    dipole(vx, vy, vz, psi, mux, muy, muz, eval, ist, par);
    mag_dipole(vx, vy, vz, psi, mx, my, mz, eval, planfw, planbw, fftwpsi, ist, par);
    
    // momentum(vx, vy, vz, psi, px_ni, py_ni, pz_ni, eval, planfw, planbw, fftwpsi, ist, par);

    rotational_strength(rs, mux, muy, muz, mx, my, mz, eval, ist);
    bethe_salpeter(bsmat, h0mat, psi, vz, mux, muy, muz, mx, my, mz, ist, par);
  
    writeSeparation(stdout);
    writeCurrentTime(stdout);

    /***********************************************************************/
    free(psi); free(potq); free(potqx); free(eval);  free(de); 
    free(ksqr); free(vx); free(vy);  free(vz); 
    free(rx); free(ry); free(rz); free(poth);  free(potl);
    free(bsmat); free(h0mat); free(mux); free(muy); free(muz);
    free(mx); free(my); free(mz); free(rs);
    free(planfw); free(planbw);

    return 0;
}

/*****************************************************************************/
