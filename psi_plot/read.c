/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/

void read_conf(double *rx, double *ry, double *rz, atm_st *atm, long ntot, FILE *pf) {
  FILE *pw;
  long i; 
  double xd, yd, zd;
  
  for (xd = yd = zd = 0.0, i = 0; i < ntot; i++) {
    fscanf(pf, "%s %lf %lf %lf", atm[i].atyp, &rx[i], &ry[i], &rz[i]);
    atm[i].natyp = assign_atom_number(atm[i].atyp);
    xd += rx[i];
    yd += ry[i];
    zd += rz[i];
  }
  xd /= (double)(ntot);
  yd /= (double)(ntot);
  zd /= (double)(ntot);
  for (i = 0; i < ntot; i++) {
    rx[i] -= xd;
    ry[i] -= yd;
    rz[i] -= zd;
  }  

  pw = fopen("conf.dat" , "w");
  fprintf(pw,"%ld\n",ntot);
  for (i = 0; i < ntot; i++) {
    // fprintf(pw, "%s %g %g %g %ld\n", atm[i].atyp, rx[i], ry[i], rz[i], atm[i].natyp);
    fprintf(pw, "%s %g %g %g\n", atm[i].atyp, rx[i], ry[i], rz[i]);
  }
  fclose(pw);
  
  return;
}

/*****************************************************************************/

void read_pot(double *vr, double *pot, long *npot, double *dr, atm_st *atm, long n, long ntype, double *a4Params) {
  FILE *pf;  
  long i, j, iscan, nloc; 
  char str[100], str2[100], atype[3];
  double *a, *b;

  if ((a = (double *) calloc(ntype, sizeof(double)))==NULL) nerror("a");
  if ((b = (double *) calloc(ntype, sizeof(double)))==NULL) nerror("b");

  a[8] = 0.64;
  a[9] = -0.384;
  a[10] = 0.64;
  a[11] = -0.20;
  b[8] = b[9] = b[10] = b[11] = 2.2287033;

  /*
      Cd = 0
      Se = 1
      In = 2
      As = 3
      Si = 4
      H  = 5
      Zn = 6
      S  = 7
      P1 = 8
      P2 = 9
      P3 = 10
      P4 = 11
      Te = 12
      Cdz = 13
      Sez = 14
      Ga = 15
      P = 16
  */

  fprintf(stdout, "ntype = %ld\n", ntype);

  for (j = 0; j < ntype*n; j++) pot[j] = vr[j] = 0;
  for (j = 0; j < ntype; j++) npot[j] = 0;
  
  for (j = 0; j < ntype; j++) {
    assign_atom_type(atype,j);
    if ((j==5) || (j==7) || (j==16)) {
        sprintf(str, "pot%c.par", atype[0]); 
        sprintf(str2, "pot%c_a4.par", atype[0]);
    }
    else if ((j <= 12) || (j==15)) {
        sprintf(str,"pot%c%c.par", atype[0], atype[1]); 
        sprintf(str2, "pot%c%c_a4.par", atype[0], atype[1]);
    }
    else {
        sprintf(str, "pot%c%c%c.par", atype[0], atype[1], atype[2]); 
        sprintf(str2, "pot%c%c%c_a4.par", atype[0], atype[1], atype[2]);
    }
    pf = fopen(str , "r");
    if (pf != NULL) {
      for (npot[j] = iscan = 0; iscan != EOF; npot[j]++) {
	    iscan = fscanf(pf, "%lg %lg", &vr[j*n+npot[j]], &pot[j*n+npot[j]]);
      }
	  fclose(pf);
      npot[j]--;
    }      
    else {
      npot[j] = npot[0];
      for (i = 0; i < npot[j]; i++) {
	    vr[j*n+i] = vr[i];
	    pot[j*n+i] = (a[j] * exp(-sqr(vr[j*n+i]) / b[j]));
      }
    }

    pf = fopen(str2, "r");
    if (pf != NULL) {
        fscanf(pf, "%lg", &a4Params[j]);
        fclose(pf);
    }
    else {
        a4Params[j] = 0.;
        printf("Warning: no %s file... setting a4 to 0.\n", str2);
    }

    /*** shift the potentials and get the r-spacing ***/
    for (i = 0; i < npot[j]; i++) {
      pot[j*n+i] -= pot[j*n+npot[j]-1];
      dr[j] = vr[j*n+1] - vr[j*n+0];
    }
      
    /*** print shifted pot ***/
    if ((j==5) || (j==7) || (j==16)) sprintf(str, "pot%c.dat", atype[0]);
    else if ((j <= 12) || (j==15)) sprintf(str, "pot%c%c.dat", atype[0], atype[1]);
    else sprintf(str, "pot%c%c%c.dat", atype[0], atype[1], atype[2]);
    pf = fopen(str , "w");
    for (i = 0; i < npot[j]; i++) fprintf (pf, "%g %g\n", vr[j*n+i], pot[j*n+i]);
    fclose(pf);
  }

  free(a); free(b);

  return;
}

/****************************************************************************/

long assign_atom_number(char atyp[3])
{
  char strerror[100];
  
  if ((atyp[0] == 'C') && (atyp[1] == 'd')  && (atyp[2] == '\0')) return(0);
  else if ((atyp[0] == 'S') && (atyp[1] == 'e') && (atyp[2] == '\0')) return(1);
  else if ((atyp[0] == 'I') && (atyp[1] == 'n') && (atyp[2] == '\0')) return(2);
  else if ((atyp[0] == 'A') && (atyp[1] == 's') && (atyp[2] == '\0')) return(3);
  else if ((atyp[0] == 'S') && (atyp[1] == 'i') && (atyp[2] == '\0')) return(4);
  else if ((atyp[0] == 'H') && (atyp[1] == '\0') && (atyp[2] == '\0')) return(5);
  else if ((atyp[0] == 'Z') && (atyp[1] == 'n') && (atyp[2] == '\0')) return(6);
  else if ((atyp[0] == 'S') && (atyp[1] == '\0') && (atyp[2] == '\0')) return(7);
  else if ((atyp[0] == 'P') && (atyp[1] == '1') && (atyp[2] == '\0')) return(8);
  else if ((atyp[0] == 'P') && (atyp[1] == '2') && (atyp[2] == '\0')) return(9);
  else if ((atyp[0] == 'P') && (atyp[1] == '3') && (atyp[2] == '\0')) return(10);
  else if ((atyp[0] == 'P') && (atyp[1] == '4') && (atyp[2] == '\0')) return(11);
  else if ((atyp[0] == 'T') && (atyp[1] == 'e') && (atyp[2] == '\0')) return(12);
  else if ((atyp[0] == 'C') && (atyp[1] == 'd') && (atyp[2] == 'z')) return(13);
  else if ((atyp[0] == 'S') && (atyp[1] == 'e') && (atyp[2] == 'z')) return(14);
  else if ((atyp[0] == 'G') && (atyp[1] == 'a') && (atyp[2] == '\0')) return(15);
  else if ((atyp[0] == 'P') && (atyp[1] == '\0') && (atyp[2] == '\0')) return(16);
  else {
    sprintf (strerror,"atom type %s not in current list",atyp);
    nerror (strerror);
  }
  return(0);
}

/****************************************************************************/

void assign_atom_type(char *atyp,long j)
{
  if (j == 0) {atyp[0] = 'C'; atyp[1] = 'd'; atyp[2] = '\0';}
  else if (j == 1) {atyp[0] = 'S'; atyp[1] = 'e'; atyp[2] = '\0';}
  else if (j == 2) {atyp[0] = 'I'; atyp[1] = 'n'; atyp[2] = '\0';}
  else if (j == 3) {atyp[0] = 'A'; atyp[1] = 's'; atyp[2] = '\0';}
  else if (j == 4) {atyp[0] = 'S'; atyp[1] = 'i'; atyp[2] = '\0';}
  else if (j == 5) {atyp[0] = 'H'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 6) {atyp[0] = 'Z'; atyp[1] = 'n'; atyp[2] = '\0';}
  else if (j == 7) {atyp[0] = 'S'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 8) {atyp[0] = 'P'; atyp[1] = '1'; atyp[2] = '\0';}
  else if (j == 9) {atyp[0] = 'P'; atyp[1] = '2'; atyp[2] = '\0';}
  else if (j == 10) {atyp[0] = 'P'; atyp[1] = '3'; atyp[2] = '\0';}
  else if (j == 11) {atyp[0] = 'P'; atyp[1] = '4'; atyp[2] = '\0';}
  else if (j == 12) {atyp[0] = 'T'; atyp[1] = 'e'; atyp[2] = '\0';}
  else if (j == 13) {atyp[0] = 'C'; atyp[1] = 'd'; atyp[2] = 'z';}
  else if (j == 14) {atyp[0] = 'S'; atyp[1] = 'e'; atyp[2] = 'z';}
  else if (j == 15) {atyp[0] = 'G'; atyp[1] = 'a'; atyp[2] = '\0';}
  else if (j == 16) {atyp[0] = 'P'; atyp[1] = '\0'; atyp[2] = '\0';}
  return;
}

/*****************************************************************************/
