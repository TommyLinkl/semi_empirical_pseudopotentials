/* This files calculates the basis vectors up to a specified KE cutoff 
    given the KE cutoff and reciprocal space lattice vectors (gVectors) */

/*****************************************************************************************/

#include "fit.h"

/*****************************************************************************************/

int calcBasisVectors(vector *tmpBV, vector *gVectors, double maxKE) {
  double minGmag = 100000.0;
  int nx, ny, nz, nbv;
  int i, j, k, idim, index;

  for (i = 0; i < NDIM; i++) {
    if (gVectors[i].mag < minGmag) minGmag = gVectors[i].mag;
  }

  nx = ny = nz = (int) (sqrt(maxKE / sqr(minGmag))) + 1;   
 
  index = 0; nbv = 0;
  for (i = -nx; i < nx + 1; i++) {
    for (j = -ny; j < ny + 1; j++) {
      for (k = -nz; k < nz + 1; k++) {
        tmpBV[index].x = i*gVectors[0].x + j*gVectors[1].x + k*gVectors[2].x;
        tmpBV[index].y = i*gVectors[0].y + j*gVectors[1].y + k*gVectors[2].y;
        tmpBV[index].z = i*gVectors[0].z + j*gVectors[1].z + k*gVectors[2].z;
        tmpBV[index].mag = sqrt(retDotProduct(tmpBV[index], tmpBV[index]));
        if (HBAR*0.5*sqr(tmpBV[index].mag)/MASS < maxKE) { 
          nbv++;
        }
        index++;
      }
    }
  }

  return nbv;
}

/*****************************************************************************************/

void finalizeBasisStates(vector *basisStates, vector *tmpBV, int nMaxBV, double maxKE) {
  FILE *pf;
  int i, nbv;

  nbv = 0;
  for (i = 0; i < nMaxBV; i++) {
    if (HBAR*0.5*sqr(tmpBV[i].mag)/MASS < maxKE) {
      basisStates[nbv].x = tmpBV[i].x;
      basisStates[nbv].y = tmpBV[i].y;
      basisStates[nbv].z = tmpBV[i].z;
      basisStates[nbv].mag = tmpBV[i].mag;
      nbv++;
    }
  }

  quickSortVectors(basisStates, 0, nbv-1);

  return;
}

/*****************************************************************************************/

int countDistinctMagnitudes(vector *vect, int nVectors) {
  int i, count;

  count = 1;
  for (i = 1; i < nVectors; i++) {
    if (vect[i].mag <= vect[i-1].mag) continue;
    else count++;
  }

  return count;
}

/*****************************************************************************************/

void quickSortVectors(vector *vectors, int l, int r) {
   int j;

   if (l < r) {
      // divide and conquer
      j = partitionVectors(vectors, l, r);
      quickSortVectors(vectors, l, j-1);
      quickSortVectors(vectors, j+1, r);
   }

   return;
}

/*****************************************************************************************/

int partitionVectors(vector *vectors, int l, int r) {
   int i, j;
   double pivot = vectors[l].mag;
   vector temp;

   i = l; j = r + 1;

   while(1) {
      do ++i; while (vectors[i].mag <= pivot && i <= r);
      do --j; while (vectors[j].mag > pivot);
      if (i >= j) break;
      temp = vectors[i]; vectors[i] = vectors[j]; vectors[j] = temp;
   }
   temp = vectors[l]; vectors[l] = vectors[j]; vectors[j] = temp;

   return j;
}

/*****************************************************************************************/