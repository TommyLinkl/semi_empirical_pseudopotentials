/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/

double interpolate(double r,double dr,double *vr,double *pot,long npot,long n,long j, double strainScaleFactor)
{
  double a, b;
  long i;

  i = (long)(r / dr);
  // printf("r, dr, n, r/dr, i = %g %g %g %g %g\n", r, dr, n, r/dr, i); 
  if (i > (n - 2)) return (0.0);

  a = strainScaleFactor * (pot[j*npot+i+1] - pot[j*npot+i]) / (vr[j*npot+i+1] - vr[j*npot+i]); 
  b = strainScaleFactor * pot[j*npot+i] - vr[j*npot+i] * a; 
  
  // printf ("a, b, r, strainScaleFactor = %g %g %g %g\n",a,b,r,strainScaleFactor);
  return(a * r + b);
}

/*****************************************************************************/
