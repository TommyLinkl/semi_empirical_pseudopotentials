#include "fd.h"

#define IM1 2147483563
#define IM2 2147483399
#define AM  (1.0 / (double)IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define REPS 1.2e-7
#define RNMX (1.0 - REPS)

/****************************************************************************/

double ran_nrc(long *idum)
{
  long    j;
  long   k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double temp;
  
  if (*idum <= 0){  /* Initialize */
    if (-(*idum) < 1) *idum = 1;    /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2 = (*idum);

    for (j = NTAB+7; j >= 0; j--){/* Load the shuffle table after 8 warm-ups */
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k * IQ1) - k * IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }

  k = (*idum) / IQ1;      /* Start Here when not initializing */
  *idum = IA1 * (*idum - k * IQ1) - k * IR1; /* Compute idum = IA1*idum % IM1*/
  if (*idum < 0) *idum += IM1;               /* without overflow by Scharge  */

  k = idum2 / IQ2;                           /* method.                      */
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2; /* Compute idum2=IA2*idum % IM2 */
  if (idum2 < 0) idum2 += IM2;               /* likewise.                    */
  
  j = iy / NDIV;          /* Will be in the range 0..NTAB-1 */
  iy = iv[j] - idum2;     /* Here idum is shuffled, idum and idum2 are */
  iv[j] = *idum;          /* combined to generate output               */

  if (iy < 1) iy += IMM1;
  if ((temp = AM * iy) > RNMX) return RNMX; /* Because users don't expect */
  else return temp;                         /* endpolong value             */
}

/****************************************************************************/

void Randomize()
{
  double tt = (double)time(0);
  long   rn = (long)tt % 1000;
  long    i;
  
  for (i = 0; i < rn; i++) random();
  return;
}

/****************************************************************************/

double ran()
{
  double     rnd;
  long static once = 0;
  
  if(!once) {
    struct timeval tv;
    struct timezone tpz;

    gettimeofday(&tv,&tpz);
    srandom(tv.tv_sec);
    once = 1;
  }
  rnd = (double)random() / 2147483647.0;

  return rnd;
}
