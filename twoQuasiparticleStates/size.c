/*****************************************************************************/
//
// This file deals with the size of a nanocrystal
// 
//
/*****************************************************************************/

#include "qp.h"

/*****************************************************************************/
// Returns the range of the double array, i.e. max(dArray[i]-dArray[j])
// returns -1.0 if nElements = 0

double retRangeOfDoubleArray(double *dArray, long nElements) {
  long i, j;
  double testRange, range = -1.0;

  for (i = 0; i < nElements; i++) {
    for (j = i+1; j < nElements; j++) {
      testRange = fabs(dArray[i]-dArray[j]);
      if (testRange > range) {
        range = testRange;
      }
    }
  }

  return range;
}

/*****************************************************************************/