/****************************************************************************/
//
//
//
/****************************************************************************/

#ifndef VECTOR_H
#define VECTOR_H

/****************************************************************************/
/* These are the library functions that are used within vector.c */

#include <math.h>

/****************************************************************************/
/* These are some common multiplication schemes */

#define MAGEPS		 1e-12
#define sqr(x)       ((x) * (x))

/****************************************************************************/
/* Structure declarations */

typedef struct vector_ {
  double x, y, z;
  double mag;
} vector;

/****************************************************************************/
/* Function declarations - public interface */

/* Functions for vector structures - in vector.c */
vector retAddedVectors(vector vect1, vector vect2);
vector retSubtractedVectors(vector vect1, vector vect2);
vector retScaledVector(vector vect, double scale);
vector retCrossProduct(vector vect1, vector vect2);
vector retNormedVector(vector vect);
vector retZeroVector(void);
double retDotProduct(vector vect1, vector vect2);
double retVectorMagnitude(vector vect);
double retDistanceBetweenPoints(vector vect1, vector vect2);
int compareVectorMagnitudes(vector vect1, vector vect2);
int countDistinctMagnitudes(vector *vect, int nVectors);
void quickSortVectors(vector *vectors, int l, int r);
int partitionVectors(vector *vectors, int l, int r);
void deepCopyVectors(vector *outputVectors, vector *inputVectors, long numVectors);
void copyDoublesToVectors(vector *vectors, double *array, long numVectors);
void copyVectorsToDoubles(double *array, vector *vectors, long numVectors);

/****************************************************************************/

#endif

/****************************************************************************/