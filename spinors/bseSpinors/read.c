#include "fd.h"

/*****************************************************************************/

void read_conf(double *rx, double *ry, double *rz, atm_st *atm, long ntot, FILE *pf)
{
  FILE *pw;
  long i; 
  double xd, yd, zd;
  
  for (xd = yd = zd = 0.0, i = 0; i < ntot; i++) {
    fscanf (pf,"%s %lf %lf %lf",atm[i].atyp,&rx[i],&ry[i],&rz[i]);
    atm[i].natyp = assign_atom_number(atm[i].atyp);
    xd += rx[i];
    yd += ry[i];
    zd += rz[i];
  }
  xd /= (double)(ntot);
  yd /= (double)(ntot);
  zd /= (double)(ntot);

  for (i = 0; i < ntot; i++){
    rx[i] -= xd;
    ry[i] -= yd;
    rz[i] -= zd;
  }

  pw = fopen("conf.dat" , "w");
  for (i = 0; i < ntot; i++) fprintf (pw,"%s %g %g %g %ld\n", atm[i].atyp, rx[i], ry[i], rz[i], atm[i].natyp);
  fclose(pw);
  
  return;
}

/*****************************************************************************/

long assign_atom_number(char atyp[3]) {
  char strerror[100];
  
  if ((atyp[0] == 'C') && (atyp[1] == 'd')  && (atyp[2] == '\0')) return(0);
  else if ((atyp[0] == 'S') && (atyp[1] == 'e') && (atyp[2] == '\0')) return(1);
  else if ((atyp[0] == 'I') && (atyp[1] == 'n') && (atyp[2] == '\0')) return(2);
  else if ((atyp[0] == 'A') && (atyp[1] == 's') && (atyp[2] == '\0')) return(3);
  else if ((atyp[0] == 'S') && (atyp[1] == 'i') && (atyp[2] == '\0')) return(4);
  else if ((atyp[0] == 'H') && (atyp[1] == '\0') && (atyp[2] == '\0'))  return(5);
  else if ((atyp[0] == 'Z') && (atyp[1] == 'n') && (atyp[2] == '\0'))  return(6);
  else if ((atyp[0] == 'S') && (atyp[1] == '\0') && (atyp[2] == '\0'))  return(7);
  else if ((atyp[0] == 'P') && (atyp[1] == '1') && (atyp[2] == '\0'))  return(8);
  else if ((atyp[0] == 'P') && (atyp[1] == '2') && (atyp[2] == '\0'))  return(9);
  else if ((atyp[0] == 'P') && (atyp[1] == '3') && (atyp[2] == '\0'))  return(10);
  else if ((atyp[0] == 'P') && (atyp[1] == '4') && (atyp[2] == '\0'))  return(11);
  else if ((atyp[0] == 'T') && (atyp[1] == 'e') && (atyp[2] == '\0')) return(12);
  else if ((atyp[0] == 'C') && (atyp[1] == 'd') && (atyp[2] == 'z')) return(13);
  else if ((atyp[0] == 'S') && (atyp[1] == 'e') && (atyp[2] == 'z')) return(14);
  else {
    sprintf (strerror,"atom type %s not in current list",atyp);
    nerror (strerror);
  }
  
  return(0);
}

/*****************************************************************************/