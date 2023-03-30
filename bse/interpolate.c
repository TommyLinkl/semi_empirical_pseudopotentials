/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/

double interpolate(double r,double dr,double *vr,double *pot,int npot,int n,int j) {
  double a, b;
  int i;

  i = (int)(r / dr);
  if (i > (n - 2)) return (0.0);

  a = (pot[j*npot+i+1] - pot[j*npot+i]) / (vr[j*npot+i+1] - vr[j*npot+i]);
  b = pot[j*npot+i] - vr[j*npot+i] * a;
  
  return(a * r + b);
}

/*****************************************************************************/