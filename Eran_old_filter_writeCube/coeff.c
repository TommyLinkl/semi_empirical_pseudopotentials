#include "fd.h"

#define filt_func(x,dt)   (sqrt(dt / PIE) * exp(-sqr(x) * dt))

void coefficient(zomplex *an,zomplex *sampoints,double *target,long ncheby,long ntarget,long nthreads,double dE,double Emin,double dt)
/*** an complex array size ncheby * ntarget. On output: an[ncheby*jtarget+jcheby] contains the newton expansion coefficients ***/ 
/*** sampoints complex array size ncheby. On output: sampoints[jcheby] contains the sampling points. ***/
/*** target double array size ntarget. Input: The target energies to filter  ***/
/*** ncheby: length of newton interpolation polynomial ***/
/*** ntarget: number of target filter states per cycle ***/
/*** nthreads: number of threads ***/
/*** dE: energy range of the Hamiltonian ***/
/*** Emin: Minimum eigenvalue of the Hamiltonian ***/
/*** dt: inverse width of the Gaussian filter ***/
{
  FILE *pf;
  double scale, res = 1.0, Smin = -2.0, Smax = 2.0, x, rho, sumre, sumim;
  long icheby, jcheby, jtarget;
  
  scale = (Smax - Smin) / dE;
  rho = samp_points(sampoints,Smin,Smax,ncheby);

  omp_set_dynamic(0);
  omp_set_num_threads(nthreads);
#pragma omp parallel for private(x,jcheby,icheby,res,sumre,sumim)
  for (jtarget = 0; jtarget < ntarget; jtarget++){
    x = (sampoints[0].re + 2.0) / scale + Emin;
    an[ncheby*jtarget+0].re = filt_func(x-target[jtarget],dt);
    an[ncheby*jtarget+0].im = 0.0;
    //fprintf (pf,"%ld %g %g\n",0,an[ncheby*jtarget+0].re,an[ncheby*jtarget+0].im);
    
    x = (sampoints[1].re + 2.0) / scale + Emin;
    an[ncheby*jtarget+1].re = (filt_func(x-target[jtarget],dt) -
		      an[ncheby*jtarget+0].re) / (sampoints[1].re - sampoints[0].re);
    an[ncheby*jtarget+1].im = (-an[ncheby*jtarget+0].im) / (sampoints[1].re - sampoints[0].re);
    //fprintf (pf,"%ld %g %g\n",1,an[ncheby*jtarget+1].re,an[ncheby*jtarget+1].im);
    for (jcheby = 2; jcheby < ncheby; jcheby++) {
      for (icheby = 1, res = 1.0, sumre = sumim = 0.0; icheby < jcheby; icheby++) {
	res *= (sampoints[jcheby].re - sampoints[icheby-1].re);
	sumre += an[ncheby*jtarget+icheby].re * res;
	sumim += an[ncheby*jtarget+icheby].im * res;
      }
      res *= (sampoints[jcheby].re - sampoints[jcheby-1].re);
      x = (sampoints[jcheby].re + 2.0) / scale + Emin;
      an[ncheby*jtarget+jcheby].re = (filt_func(x-target[jtarget],dt) -
			an[ncheby*jtarget+0].re - sumre) / res;
      an[ncheby*jtarget+jcheby].im = (-an[ncheby*jtarget+0].im - sumim) / res;
      //fprintf (pf,"%ld %g %g\n",jcheby,an[ncheby*jtarget+jcheby].re,an[ncheby*jtarget+jcheby].im);
    }
  }
  return;
}

/**************************************************************************/

double samp_points(zomplex *point,double min,double max,long nc)
{
  long j,k,m, nc3=32*nc, jrnd, kmax, imfrac;
  double dsre, dsim, rnd, fkmax, fk, minim, maxim, range = 0.0, *veca, del;
  zomplex *samp;

  samp = (zomplex*)calloc(nc3,sizeof(zomplex));  
  veca = (double*)calloc(nc3,sizeof(double));  

  minim = range * min;
  maxim = range * max;
  imfrac = nc3 / 8;
  dsre = (max-min)/(double)(nc3-1-imfrac);
  dsim = (maxim-minim)/(double)(imfrac-1);
  for (j = 0; j < nc3; j++) samp[j].re = samp[j].im = 0.0;
  for (j = 0; j < nc3-imfrac; j++) samp[j].re = min + (double)(j) * dsre;
  for (j = nc3-imfrac; j < nc3; j++) samp[j].im = minim + (double)(j-nc3+imfrac) * dsim;

  jrnd = 0;
  point[0] = samp[jrnd];

  for (k = 0; k < nc3; k++){
    del = sqr(samp[k].re-point[0].re) + sqr(samp[k].im-point[0].im);
    if (del < 1.0e-10) veca[k] = -1.0e30;
    else veca[k] = log(del);
  }
    
  for (j = 1; j < nc; j++){
    kmax = 0; fkmax = -2.0e30;
    for (k = 0; k < nc3; k++){
      if (veca[k] > fkmax){
	fkmax = veca[k];
	kmax = k;
      }
    }
    point[j] = samp[kmax];
    for (k = 0; k < nc3; k++){
      del = sqr(samp[k].re-point[j].re) + sqr(samp[k].im-point[j].im);
      if (del < 1.0e-10) veca[k] = -1.0e30;
      else veca[k] += log(del);
    }
    
  }
  for (fk = 1.0, j = 0; j < nc; j++) fk *= (sqr(point[j].re) + sqr(point[j].im));
  fk = sqrt(fk);
  fk = pow(fk,1.0/(double)(nc));
  /*for (j = 0; j < nc; j++) {
    point[j].re /= fk;
    point[j].im /= fk;
    }*/
  free(samp); free(veca);
  return (fk);
}

/**************************************************************************/

