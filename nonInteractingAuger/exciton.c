#include "fd.h"

/*****************************************************************************/

double exciton_energy(double *Cbs,double *Hbs,lng_st ist)
{
  double sum; int i, j;

  for (sum = 0.0, i = 0; i < ist.msbs2; i++)
    for (j = 0; j < ist.msbs2; j++)
      sum += Cbs[j] * Cbs[i] * Hbs[i*ist.msbs2+j];
  
  return(sum);
}

/*****************************************************************************/