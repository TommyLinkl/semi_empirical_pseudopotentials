/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

/****************************************************************************/

#define filt_func(x,dt)   (sqrt(dt / PIE) * exp(-sqr(x) * dt))
//#define filt_func(x,dt)   (0.04 / (PIE * (sqr(x) + sqr(0.04))))

/****************************************************************************/

void coefficient(zomplex *an,double *samp,double *El,par_st par,long_st ist)
{
  FILE *pf;
  zomplex *samploc;
  double scale, res = 1.0, Smin = -2.0, Smax = 2.0, x, rho, sumre, sumim;
  long i, j, ie;
  
  scale = (Smax - Smin) / par.dE;
  //chebyshev_reordered(samp,Smin,Smax,ist.nc);
  samploc = (zomplex*)calloc(ist.nc,sizeof(zomplex));
  rho = samp_points_ashkenazy(samploc,Smin,Smax,ist.nc);
  for (j = 0; j < ist.nc; j++) samp[j] = samploc[j].re;

  omp_set_dynamic(0);
  omp_set_num_threads(ist.nthreads);
#pragma omp parallel for private(x,j,i,res,sumre,sumim)
  for (ie = 0; ie < ist.ms; ie++){
    x = (samp[0]+2.0)/scale+par.Vmin;
    an[ist.nc*ie+0].re = filt_func(x-El[ie],par.dt);
    an[ist.nc*ie+0].im = 0.0;
    //fprintf (pf,"%ld %g %g\n",0,an[ist.nc*ie+0].re,an[ist.nc*ie+0].im);
    
    x = (samp[1]+2.0)/scale+par.Vmin;
    an[ist.nc*ie+1].re = (filt_func(x-El[ie],par.dt) -
		      an[ist.nc*ie+0].re) / (samp[1] - samp[0]);
    an[ist.nc*ie+1].im = (-an[ist.nc*ie+0].im) / (samp[1] - samp[0]);
    //fprintf (pf,"%ld %g %g\n",1,an[ist.nc*ie+1].re,an[ist.nc*ie+1].im);
    for (j = 2; j < ist.nc; j++) {
      for (i = 1, res = 1.0, sumre = sumim = 0.0; i < j; i++) {
	res *= (samp[j] - samp[i-1]);
	sumre += an[ist.nc*ie+i].re * res;
	sumim += an[ist.nc*ie+i].im * res;
      }
      res *= (samp[j] - samp[j-1]);
      x = (samp[j]+2.0)/scale+par.Vmin;
      an[ist.nc*ie+j].re = (filt_func(x-El[ie],par.dt) -
			an[ist.nc*ie+0].re - sumre) / res;
      an[ist.nc*ie+j].im = (-an[ist.nc*ie+0].im - sumim) / res;
      //fprintf (pf,"%ld %g %g\n",j,an[ist.nc*ie+j].re,an[ist.nc*ie+j].im);
    }
  }
  pf = fopen("coeff.dat" , "w");
  fclose(pf);

  check_function(an,samploc,ist,par,El[0]);
  free(samploc);
  return;
}

/****************************************************************************/

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

/****************************************************************************/

double samp_points_ashkenazy(zomplex *point,double min,double max,long nc)
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

/****************************************************************************/

void check_function(zomplex *an,zomplex *samp,long_st ist,par_st par,double El)
{
  FILE *pf; long j, i;
  zomplex f, x, xn, xm1, ctmp;
  double dx = 0.01, xunsc;

  pf = fopen("func.dat" , "w");
  for (x.im = 0.0, x.re = -2.0; x.re < 2.0; x.re += dx){
    xn.re = 1.0;
    xn.im = 0.0;
    for (f = an[0], j = 1; j < ist.nc; j++){
      xm1.re = x.re - samp[j-1].re;
      xm1.im = x.im - samp[j-1].im;
      cmul(xm1,xn,xn);
      cmul(an[j],xn,ctmp);
      
      f.re += ctmp.re;
      f.im += ctmp.im;
    }
    xunsc = ((x.re+2.0)* par.dE / 4.0 + par.Vmin);

    fprintf (pf,"%g %g %g\n",xunsc,f.re,filt_func(xunsc-El,par.dt));
  }
  fclose(pf);
  return;
}

/****************************************************************************/