#include "ar.h"

/*****************************************************************************/

#define filt_func(x,dt)   (sqrt(dt / PIE) * exp(-sqr(x) * dt))

/*****************************************************************************/

void coefficient(zomplex *an, double *samp, long nc,
		 long M, double dt, double dE, double Emin, double *El)
{
  FILE *pf;
  zomplex sum, *samploc;
  double scale, res = 1.0, Smin = -2.0, Smax = 2.0, x, rho;
  long i, j, ie;
  
  scale = (Smax - Smin) / dE;
  //chebyshev_reordered(samp,Smin,Smax,nc);
  if ((samploc = (zomplex*)calloc(nc,sizeof(zomplex)))==NULL)nerror("smaploc");;
  rho = samp_points_ashkenazy(samploc,Smin,Smax,nc);
  for (j = 0; j < nc; j++) samp[j] = samploc[j].re;
  free(samploc);

  for (pf = fopen("coeff.dat" , "w"), ie = 0; ie < M; ie++){
    x = (samp[0]+2.0)/scale+Emin;
    an[nc*ie+0].re = filt_func(x-El[ie],dt);
    an[nc*ie+0].im = 0.0;
    fprintf (pf,"%ld %g %g\n",0,an[nc*ie+0].re,an[nc*ie+0].im);
    
    x = (samp[1]+2.0)/scale+Emin;
    an[nc*ie+1].re = (filt_func(x-El[ie],dt) - 
		      an[nc*ie+0].re) / (samp[1] - samp[0]);
    an[nc*ie+1].im = (-an[nc*ie+0].im) / (samp[1] - samp[0]);
    fprintf (pf,"%ld %g %g\n",1,an[nc*ie+1].re,an[nc*ie+1].im);
    for (j = 2; j < nc; j++) {
      for (i = 1, res = 1.0, sum.re = sum.im = 0.0; i < j; i++) {
	res *= (samp[j] - samp[i-1]);
	sum.re += an[nc*ie+i].re * res;
	sum.im += an[nc*ie+i].im * res;
      }
      res *= (samp[j] - samp[j-1]);
      x = (samp[j]+2.0)/scale+Emin;
      an[nc*ie+j].re = (filt_func(x-El[ie],dt) - 
			an[nc*ie+0].re - sum.re) / res;
      an[nc*ie+j].im = (-an[nc*ie+0].im - sum.im) / res;
      fprintf (pf,"%ld %g %g\n",j,an[nc*ie+j].re,an[nc*ie+j].im);
    }
  }
  fclose(pf);
  return;
}

/*****************************************************************************/

void chebyshev_reordered(double *point,double min,double max,long nc)
{
  long j,k,m;

  point[0] = 1;
  for (k = 1, m = nc/2; k < nc; k *= 2, m /= 2)
    for (j = 0; j < nc - k; j++)
      point[k+j] = point[j] + m;
  for (j = 0; j < nc; j++)
    point[j]= (max+min)/2.0 + (max-min)/2.0*cos(PIE*(point[j]-0.5)/(nc));
  return;
}

/*****************************************************************************/

double samp_points_ashkenazy(zomplex *point,double min,double max,long nc)
{
  long j,k,m, nc3=32*nc, jrnd, kmax, imfrac;
  double dsre, dsim, rnd, fkmax, fk, minim, maxim, range = 0.0, *veca, del;
  zomplex *samp;

  if ((samp = (zomplex*)calloc(nc3,sizeof(zomplex)))==NULL)nerror("samp"); 
  if ((veca = (double*)calloc(nc3,sizeof(double)))==NULL)nerror("veca");

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

/*****************************************************************************/
