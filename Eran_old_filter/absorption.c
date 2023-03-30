#include "fd.h"

/****************************************************************************/

void absorption_spectrum(double *vx, double *vy, double *vz, double *psi, double *eval, double *sige, long_st ist, par_st par) {
  FILE *pf1, *pf2, *pf3, *pf; 
  long a, b, i, jx, jy, jz, jgrid, jyz; 
  long num_holes, num_electrons, num_excitons, hole_index, electron_index, *eigenstate_index_list;
  double sumx, sumy, sumz;
  double dz, dx, dy, ev, os, tmp;
  double *mux, *muy, *muz;

  /*** get the number of homo and lumo states ***/  
  num_holes = num_eigenstates_energy_range(eval, sige, -100.0, -0.18, ist.mstot);
  num_electrons = num_eigenstates_energy_range(eval, sige, -0.18, 100.0, ist.mstot);
  num_excitons = num_holes * num_electrons;

  /*** allocate memory for list of the indexes of all the eigenstates ***/
  eigenstate_index_list = (long *)calloc(num_holes+num_electrons, sizeof(long));
  printf("Number of hole eigenstates = %ld\n", num_holes);
  printf("Number of electron eigenstates = %ld\n", num_electrons);
  printf("Number of excitonic states = %ld\n", num_excitons);
  
  /*** fill up eigenstate_index_list with the index of the eigenstates ***/
  get_eigenstate_list(eigenstate_index_list, eval, sige, -100.0, 100.0, ist.mstot);

  /*** allocate memory for mux, muy, and muz ***/
  mux = (double *)calloc(num_excitons, sizeof(double));
  muy = (double *)calloc(num_excitons, sizeof(double));
  muz = (double *)calloc(num_excitons, sizeof(double));

  /*** get file pointers to print oscillator strength files ***/
  pf = fopen("OS0.dat", "w"); 
  pf1 = fopen("ux.dat", "w"); 
  pf2 = fopen("uy.dat", "w"); 
  pf3 = fopen("uz.dat", "w");
  
  for (i = 0; i < num_holes; i++, fprintf(pf1,"\n"),fprintf(pf2,"\n"),fprintf(pf3,"\n")) {
    for (a = 0; a < num_electrons; a++) {
      hole_index = eigenstate_index_list[i];
      electron_index = eigenstate_index_list[num_holes+a];
      for (sumx = sumy = sumz = 0.0, jz = 0; jz < ist.nz; jz++) {
        dz = vz[jz];
        for (jy = 0; jy < ist.ny; jy++) {
          dy = vy[jy];
          jyz = ist.nx * (ist.ny * jz + jy);
          for (jx = 0; jx < ist.nx; jx++) {
            dx = vx[jx];
            jgrid = jyz + jx;
       	    tmp = psi[hole_index*ist.ngrid+jgrid] * psi[electron_index*ist.ngrid+jgrid] * par.dv;
       	    sumx +=  dx*tmp;
       	    sumy +=  dy*tmp;
       	    sumz +=  dz*tmp;
          }
        }
      }
      mux[i*num_electrons+a] =  sumx;
      muy[i*num_electrons+a] =  sumy;
      muz[i*num_electrons+a] =  sumz;
      fprintf(pf1, "%ld %ld %g\n", i, a, sumx);
      fprintf(pf2, "%ld %ld %g\n", i, a, sumy);
      fprintf(pf3, "%ld %ld %g\n", i, a, sumz);

      os=(sqr(mux[i*num_electrons+a]) + sqr(muy[i*num_electrons+a]) + sqr(muz[i*num_electrons+a]));
      ev = (eval[electron_index] - eval[hole_index]) * AUTOEV;
      fprintf(pf, "%.10f %.8f %.8f %.10f %.10f %.10f\n", sqrt(os), ev, 
        AUTONS*((4.0/3.0)*(ev*ev*ev)*os/(137.*137.*137.)),
	      sqr(mux[i*num_electrons+a]),
	      sqr(muy[i*num_electrons+a]),
	      sqr(muz[i*num_electrons+a])
	  );
    }
  }
  fclose(pf); fclose(pf1); fclose(pf2); fclose(pf3);

  // Free dynamically allocated memory
  free(mux); free(muy); free(muz); 
  free(eigenstate_index_list);

  return;
}

/****************************************************************************/

void fill_direction_averaged_array(zomplex *ave_array) {
  zomplex xx, xy, xz, yy, yz, zz; // the 6 independent entries
  
  ave_array[0].re = xx.re = 1.0;
  ave_array[0].im = xx.im = 1.0;
  ave_array[1].re = ave_array[3].re = xy.re = 1.0;
  ave_array[1].im = xy.im = 1.0;
  ave_array[3].im = -1.0*xy.im;
  ave_array[2].re = ave_array[6].re = xz.re = 1.0;
  ave_array[2].im = xz.im = 1.0;
  ave_array[6].im = -1.0*xz.im; 
  ave_array[4].re = yy.re = 1.0;
  ave_array[4].im = yy.im = 1.0;
  ave_array[5].re = ave_array[7].re = yz.re = 1.0;
  ave_array[5].im = yz.im = 1.0;
  ave_array[7].im = -1.0*yz.im;
  ave_array[8].re = zz.re = 1.0;
  ave_array[8].im = zz.im = 1.0;

  return;
}

/****************************************************************************/
