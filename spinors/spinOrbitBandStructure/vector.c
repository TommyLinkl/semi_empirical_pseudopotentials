/****************************************************************************/
/* This file contains the transformations for the vector structure */

#include "fit.h"

/****************************************************************************/

vector retAddedVectors(vector vect1, vector vect2) {
  vect1.x += vect2.x; vect1.y += vect2.y; vect1.z += vect2.z;
  vect1.mag = sqrt(retDotProduct(vect1, vect1));

  return vect1;
}

/****************************************************************************/

vector retSubtractedVectors(vector vect1, vector vect2) {
  vect1.x -= vect2.x; vect1.y -= vect2.y; vect1.z -= vect2.z;
  vect1.mag = sqrt(retDotProduct(vect1, vect1));
 
  return vect1;
}

/****************************************************************************/

vector retScaledVector(vector vect, double scale) {
  vect.x *= scale; vect.y *= scale; vect.z *= scale;
  vect.mag = sqrt(retDotProduct(vect, vect));

  return vect;
}

/****************************************************************************/

vector retCrossProduct(vector vect1, vector vect2) {
  vector crossProduct;

  crossProduct.x = vect1.y * vect2.z - vect1.z * vect2.y;
  crossProduct.y = vect1.z * vect2.x - vect1.x * vect2.z;
  crossProduct.z = vect1.x * vect2.y - vect1.y * vect2.x;
  crossProduct.mag = sqrt(retDotProduct(crossProduct, crossProduct));

  return crossProduct;
}

/****************************************************************************/

double retDotProduct(vector vect1, vector vect2) {
  return (vect1.x * vect2.x + vect1.y * vect2.y + vect1.z * vect2.z);
}

/****************************************************************************/

int compareVectorMagnitudes(vector vect1, vector vect2) {
  if (vect1.mag <= vect2.mag) return 1;
  else if (vect1.mag > vect2.mag) return 0;
  else return -1;
}

/****************************************************************************/

vector retZeroVector() {
  vector zeroVect;
  zeroVect.x = 0.0; zeroVect.y = 0.0; zeroVect.z = 0.0; zeroVect.mag = 0.0;
  return zeroVect;
}

/****************************************************************************/
