/****************************************************************************/

#include "fit.h"

/****************************************************************************/
/*	THIS PROGRAM RETURNS A NUMBER ("rnd") 				    */
/*	BETWEEN -1.0 and 1.0 .						    */
/****************************************************************************/

double randNumber() {
  double rnd;
  static int once = 0;
  
  if(!once) {
    struct timeval tv;
    gettimeofday(&tv, 0);
    srandom(tv.tv_sec);
    once = 1;
  }
  rnd = (double)random() / (double)2147483647;
  
  if (rnd > (double)1) {
    fprintf (stderr, "myrandom rnd = %f", rnd);
    exit(1);
  }
  
  if (rnd < (double)0) {
    fprintf (stderr, "myrandom rnd = %f", rnd);
    exit(1);
  }

  return (2.0 * rnd - 1.0);
}

/****************************************************************************/
