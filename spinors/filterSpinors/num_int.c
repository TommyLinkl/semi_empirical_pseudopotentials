#include "fd.h"

/*****************************************************************************/
// solves the differential equation of the form:
// u"(r) + p(r)u'(r)+ q(r)u(r) + s(r) = 0
// and stores the result in u

void gen_numerov_solver(double *u, double *p, double *q, double *s, double *r, double u0, double u1, int istart, int iend) {
	FILE *pf;
	int n, numSteps = abs(iend-istart);
	double aPrime, a2prime, h_1, h_2, h = r[1]-r[0]; // the step size
	double *a, *z, *P, *F, *A, *B; 
	double a0, z0, z1;
	
	a = (double *) calloc(numSteps, sizeof(double));	// a'(r) = P(r)*F(a(r)), runge-kutta
	z = (double *) calloc(numSteps, sizeof(double));	// z"(r) + A(r)z(r) + B(r) = 0, numerov
	P = (double *) calloc(numSteps, sizeof(double));	// P(r) = -0.5*p(r) - used in runge-kutta (dirac eqn specific)
	F = (double *) calloc(numSteps, sizeof(double));	// F(a(r)) = a(r) - used in runge-kutta (dirac eqn specific)
	A = (double *) calloc(numSteps, sizeof(double));	// A(r) = (a"(r)+p(r)a'(r)+q(r)a(r))/a(r)
	B = (double *) calloc(numSteps, sizeof(double));	// B(r) = s(r)/a(r)

	h_1 = 1.0/h;
	h_2 = 1.0/(h*h);	

	// calculate P(r) then call runge-kutta solver to get a(r)
	calc_P_rk_dirac(P, p, istart, iend);
	a0 = 0.01;
	runge_kutta_solver(a, P, F, r, a0, h, istart, iend);

	// calculate A(r) and B(r) then call numerov solver to get z(r)
	for (n = 0; n < numSteps; n++) {
		if (n > 0 && n < numSteps-1) aPrime = (a[n+1]-a[n-1])*0.5*h_1;
		else if (! n) aPrime = (a[n+1]-a[n])*h_1; 
		else aPrime = (a[n]-a[n-1])*h_1;
		if (n > 0 && n < numSteps-1) a2prime = (a[n+1]-2.0*a[n]+a[n-1])*h_2;
		else if (n == 0) a2prime = (a[n+2]-2.0*a[n+1]+a[n])*h_2;
		else a2prime = (a[n]-2.0*a[n-1]+a[n-2])*h_2;
		if (istart < iend) {
			A[n] = (a2prime+p[istart+n]*aPrime+q[istart+n]*a[n])/(a[n]+EPSDX);
			B[n] = s[istart+n]/(a[n]+EPSDX);
		}
		else {
			A[n] = (a2prime+p[iend+n]*aPrime+q[iend+n]*a[n])/(a[n]+EPSDX);
			B[n] = s[iend+n]/(a[n]+EPSDX);
		}
	}
	z0 = u0/a0;
	z1 = u1/a0;
	numerov_solver(z, A, B, z0, z1, h, istart, iend);

	// calculate u(r)=a(r)*z(r) and store it in the inputed double pointer/ array u 
	if (istart < iend) for (n = 0; n < numSteps; n++) u[istart+n] = a[n]*z[n];	
	else for (n = 0; n < numSteps; n++) u[iend+n] = a[n]*z[n];	

	if (istart < iend) {
		pf = fopen("NumInt.dat", "w");
		for (n = 0; n < numSteps; n++) fprintf(pf, "%.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", q[n], p[n], P[n], A[n], a[n], z[n], u[n]); 
		fclose(pf);
	}

	free(P); free(F);
	free(A); free(B);
	free(a); free(z);

	return;
}

/*****************************************************************************/
// solves the differential equation of the form:
// y"(r) + k(r)y(r) + b(r) = 0 
// and stores the result in y(r)

void numerov_solver(double *y, double *k, double *b, double y0, double y1, double h, int istart, int iend) {
	// TODO: Add in b(r) - not needed in the relativistic Schrodinger equation solver though
	int n, numSteps = abs(iend-istart);
	double one_twelth_h2 = h*h/12.0;

	if (istart < iend) {
		y[0] = y0;
		y[1] = y1;
		for (n = 1; n < numSteps-1; n++) {
			y[n+1] = 2.0*(1.0-5*one_twelth_h2*k[n])*y[n]-(1+one_twelth_h2*k[n-1])*y[n-1];
			y[n+1] /= (1+one_twelth_h2*k[n+1]);
		}
	}
	else {
		y[numSteps-1] = y0;
		y[numSteps-2] = y1;
		for (n = numSteps-2; n > 0; n--) {
                        y[n-1] = 2.0*(1.0-5*one_twelth_h2*k[n])*y[n]-(1+one_twelth_h2*k[n+1])*y[n+1];
                        y[n-1] /= (1+one_twelth_h2*k[n-1]);
                }
	}

	return;
}

/*****************************************************************************/
// solves the differential equation of the form:
// a'(r) = p(r)*f(a(r)) = p(r)*a(r) - last expression for relativistc Schrodinger eqn
// and stores the result in a(r)

void runge_kutta_solver(double *a, double *p, double *f, double *r, double a0, double h, int istart, int iend) {
	int n, numSteps = abs(iend-istart);
	double k1, k2, k3, k4;
	double one_sixth = 1.0/6.0;

	if (istart < iend) {
		a[0] = a0;	// initial condition
		for (n = 0; n < numSteps-1; n++) {
			k1 = h * p[n] * a[n];
			k2 = h * calc_average(p, n, n+1) * (a[n] + k1*0.5);
			k3 = h * calc_average(p, n, n+1) * (a[n] + k2*0.5);
			k4 = h * p[n+1] * (a[n] + k3);
			a[n+1] = a[n] + one_sixth*(k1 + 2.0*k2 + 2.0*k3 + k4);
		}
	}
	else {
		a[numSteps-1] = a0;	// initial condition
		for (n = numSteps-1; n > 0; n--) {
			k1 = h * p[n] * a[n];
			k2 = h * calc_average(p, n-1, n) * (a[n] + k1*0.5);
			k3 = h * calc_average(p, n-1, n) * (a[n] + k2*0.5);
			k4 = h * p[n-1] * (a[n] + k3);
			a[n-1] = a[n] + one_sixth*(k1 + 2.0*k2 + 2.0*k3 + k4);
		}
	}

	return;
}

/*****************************************************************************/
// this routine is specific to the relativistic Schrodinger equation
// P[n] = -0.5*p[n], P[n] sent to runge_kutta_solver

void calc_P_rk_dirac(double *P, double *p, int istart, int iend) {
	int n, tmp, numSteps = abs(iend-istart);
	if (istart < iend) tmp = istart;
	else tmp = iend;

	for (n = 0; n < numSteps; n++) P[n] = -0.5*p[tmp+n];

	return;
}

/*****************************************************************************/
// simply calculates the average for the slice of the array 
// starting at istart and ending at iend (including both endpoints)

double calc_average(double *a, int istart, int iend) {
	int i;
	double ave = 0.0; double norm = 0.0;

	for (i = istart; i < iend+1; i++) {
		ave += a[i];
		norm += 1.0;	
	}

	return (ave/norm);
}

/*****************************************************************************/

long count_nodes(double *a, long len) {
  long i, numNodes = 0;
  
  for (i = 1; i < len; i++) if (a[i-1]*a[i] <  0.0) numNodes++;
  
  return numNodes;
}

/*****************************************************************************/

double min_element(double *a, long len) {
  int i; double min = a[0];
  
  for (i = 0; i < len; i++) if ( a[i] < min ) min = a[i];
  
  return min;
}

/*****************************************************************************/
