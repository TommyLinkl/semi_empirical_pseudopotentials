#include "ar.h"

/*****************************************************************************/

double interpolate(double r,double dr,double *vr,double *pot,long nPot,long n,long j)
{
  double a, b;
  long i;

  i = (long)(r / dr);
  if (i > (n - 2)) return (0.0);

  a = (pot[j*nPot+i+1] - pot[j*nPot+i]) / (vr[j*nPot+i+1] - vr[j*nPot+i]);
  b = pot[j*nPot+i] - vr[j*nPot+i] * a;
  
  return(a * r + b);
}

/*****************************************************************************/