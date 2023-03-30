/*****************************************************************************/

#include "qp.h"

/*****************************************************************************/
/* gets largest distance between (x, y, or z) coordinates of atoms in system */
/*****************************************************************************/
double get_dot_ligand_size_z(double *rz, long n) {
  long i, j;
  double dr2, dis2;
  
  for (dr2 = 0.0, i = 0; i < n-1; i++) {
    for (j = i+1; j < n; j++) {
      dis2 = sqr(rz[i]-rz[j]);
      if (dis2 > dr2) {
        dr2 = dis2;
      }
    }
  }
  return(sqrt(dr2));
}

/*****************************************************************************/
