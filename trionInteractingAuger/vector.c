/****************************************************************************/

#include "vector.h"

/****************************************************************************/
vector retAddedVectors(vector vect1, vector vect2) {
  vect1.x += vect2.x; vect1.y += vect2.y; vect1.z += vect2.z;
  vect1.mag = retVectorMagnitude(vect1);

  return vect1;
}

/****************************************************************************/
vector retSubtractedVectors(vector vect1, vector vect2) {
  vect1.x -= vect2.x; vect1.y -= vect2.y; vect1.z -= vect2.z;
  vect1.mag = retVectorMagnitude(vect1);
 
  return vect1;
}

/****************************************************************************/
vector retScaledVector(vector vect, double scale) {
  vect.x *= scale; vect.y *= scale; vect.z *= scale;
  vect.mag = retVectorMagnitude(vect);

  return vect;
}

/****************************************************************************/
vector retElementWiseVectorMultiplication(vector vect1, vector vect2) {
  vector elementWiseMultVector;

  elementWiseMultVector.x = vect1.x * vect2.x;
  elementWiseMultVector.y = vect1.y * vect2.y;
  elementWiseMultVector.z = vect1.z * vect2.z;
  elementWiseMultVector.mag = retVectorMagnitude(elementWiseMultVector);

  return elementWiseMultVector;
}

/****************************************************************************/
vector retCrossProduct(vector vect1, vector vect2) {
  vector crossProduct;

  crossProduct.x = vect1.y * vect2.z - vect1.z * vect2.y;
  crossProduct.y = vect1.z * vect2.x - vect1.x * vect2.z;
  crossProduct.z = vect1.x * vect2.y - vect1.y * vect2.x;
  crossProduct.mag = retVectorMagnitude(crossProduct);

  return crossProduct;
}

/****************************************************************************/
vector retNormedVector(vector vect) {
  if (vect.mag > MAGEPS) return retScaledVector(vect, 1.0/vect.mag);
  else return retZeroVector();
}

/****************************************************************************/
vector retZeroVector() {
  vector zeroVect;
  zeroVect.x = 0.0; zeroVect.y = 0.0; zeroVect.z = 0.0; zeroVect.mag = 0.0;

  return zeroVect;
}

/****************************************************************************/
double retDotProduct(vector vect1, vector vect2) {
  return (vect1.x * vect2.x + vect1.y * vect2.y + vect1.z * vect2.z);
}

/****************************************************************************/
double retVectorMagnitude(vector vect) {
  return sqrt(retDotProduct(vect, vect));
}

/****************************************************************************/
double retDistanceBetweenPoints(vector vect1, vector vect2) {
  vector diffVect;
  diffVect = retSubtractedVectors(vect1, vect2);

  return diffVect.mag;
}

/****************************************************************************/
int compareVectorMagnitudes(vector vect1, vector vect2) {
  if (vect1.mag <= vect2.mag) return 1;
  else if (vect1.mag > vect2.mag) return 0;
  else return -1;
}

/****************************************************************************/
// Current implementation assumes that the vectors are sorted by their magnitudes 
// already (i.e., quickSortVectors has been called prior to this function)
/****************************************************************************/
int countDistinctMagnitudes(vector *vect, int nVectors) {
  int i, count;

  count = 1;
  for (i = 1; i < nVectors; i++) {
    if (vect[i].mag <= vect[i-1].mag) continue;
    else count++;
  }

  return count;
}

/****************************************************************************/
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

/****************************************************************************/
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

/****************************************************************************/
void deepCopyVectors(vector *outputVectors, vector *inputVectors, long numVectors) {
  long i;

  for (i = 0; i < numVectors; i++) {
    outputVectors[i].x = inputVectors[i].x;
    outputVectors[i].y = inputVectors[i].y;
    outputVectors[i].z = inputVectors[i].z;
    outputVectors[i].mag = retVectorMagnitude(inputVectors[i]);
  }

  return;
}

/****************************************************************************/
void copyDoublesToVectors(vector *vectors, double *array, long numVectors) {
  long i;

  for (i = 0; i < numVectors; i++) {
    vectors[i].x = array[3*i];
    vectors[i].y = array[3*i+1];
    vectors[i].z = array[3*i+2];
    vectors[i].mag = retVectorMagnitude(vectors[i]);
  }

  return;
}

/****************************************************************************/
void copyVectorsToDoubles(double *array, vector *vectors, long numVectors) {
  long i;

  for (i = 0; i < numVectors; i++) {
    array[3*i] = vectors[i].x;
    array[3*i+1] = vectors[i].y;
    array[3*i+2] = vectors[i].z;
  }

  return;
}

/****************************************************************************/
