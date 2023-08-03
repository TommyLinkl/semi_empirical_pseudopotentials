#include "fd.h"


void slice_psi(long_st ist, par_st par, int n, int m)
{
  double *psi, *sige;
  double tmpd;
  long tmpi,i,sig_index;
  FILE *psig, *ppsi, *ppsitot;
  psi = (double*)calloc(m*ist.ngrid,sizeof(double));
  /*** the filtered energies ***/
  if ((sige = (double*)calloc(ist.mstot,sizeof(double)))==NULL)nerror("sige");
  if(access("eval.dat", F_OK) != -1 ){
     psig = fopen("eval.dat","r");
     for (i = 0; i < ist.mstot; i++) {fscanf(psig, "%ld %lg %lg", &tmpi, &tmpd, &sige[i]);}
      fclose(psig);
  }
  else {
    printf("\n\nNo eval.dat file detected in current working directory - the program is exiting!!!\n\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
  }

  if(access("psi.dat", F_OK) != -1 ){
    ppsi = fopen("psi.dat","r");
    if (sige[n]>par.sigma_e) {printf("This is not a valid initial states to read in with sig=%lg, change ist.n!!! \n\n",sige[n]); exit(EXIT_FAILURE);}
    fseek(ppsi,(long) ist.ngrid*n*sizeof(double),SEEK_SET);
    printf ("reading the first valid states %d from psi.dat\n",n);

    sig_index = 0;
    for (i=0; i < m; i++){
      while (sige[n+sig_index]>par.sigma_e){
        printf("inside while loop, i=%ld, sig_index=%ld\n",i,sig_index);
        /*** for invalid states ***/
        fseek(ppsi,(long) ist.ngrid*sizeof(double),SEEK_CUR); // Move the pointer of ppsi
        sig_index ++;
        if (sig_index>100) {printf("Check the input, too many invalid states...\n"); exit(EXIT_FAILURE);}
      }
      fread (&psi[i*ist.ngrid],sizeof(double),ist.ngrid,ppsi);
      sig_index ++;
    }

    printf ("done reading total %d valid states\n",m);
    fclose(ppsi);
  }
  else {
    printf("\n\nNo psi.dat file detected in current working directory - the program is exiting!!!\n\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
  }


  ppsitot = fopen("psitot.dat" , "w+");
  fwrite (psi,sizeof(double),m*ist.ngrid,ppsitot);
  fclose(ppsitot);
  printf("Slicing completed. Saved psitot file!\n");

  /*** Free memeory ***/
  free(psi);

}













