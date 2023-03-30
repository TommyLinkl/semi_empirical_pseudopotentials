/*****************************************************************************/
// File contains functions that involve calculating the total size of the system 

#include "fd.h"

/*****************************************************************************/
// returns the maximum distance between 2 atoms in the nanocrystal where their
// positions are definted by (rx,ry,rz) and there are n of them

double get_dot_ligand_size(double *rx, double *ry, double *rz, long n) 
{
  long i, j;
  double dr2, dis2;
  
  for (dr2 = 0.0, i = 0; i < n-1; i++) for (j = i+1; j < n; j++) {
    dis2 = sqr(rx[i]-rx[j]) + sqr(ry[i]-ry[j]) + sqr(rz[i]-rz[j]);
    if (dis2 > dr2) dr2 = dis2;
  }

  return(sqrt(dr2));
}

/***************************************************************************/
// returns the maximum absolute value of the double z 
// within the first n doubles in rz

double get_dot_ligand_size_z(double *rz, long n)
{
  long i, j;
  double dr2, dis2;
  
  for (dr2 = 0.0, i = 0; i < n-1; i++) for (j = i+1; j < n; j++) {
      dis2 = sqr(rz[i]-rz[j]);
      if (dis2 > dr2) dr2 = dis2;
    }

  return(sqrt(dr2));
}

/*****************************************************************************/