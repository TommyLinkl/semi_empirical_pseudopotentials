/*****************************************************************************************/

/* This file calculates and prints the real-space local potential for each atom */

/*****************************************************************************************/

#include "fit.h"

/*****************************************************************************************/

#define NGRID 1024
#define QMAX  5.0
#define QMIN  0

/*****************************************************************************************/
// This function calculates and then prints the real-space potential for each atom
// by Fourier transforming the q-space potentials on a discrete uniform grid
// rSpacePot(r) = 4*pi/(2pi)^3 * sum(qSpacePot(q) * sin(qr)/(qr) * q^2 * dq) 
// where the sin(qr)/q is the spherical Bessel function for l=0

void calcRealSpacePotential(atom *atoms, param params) {
	FILE *pf;
	int i, j, k, ir;
	double *qSpacePot, *rSpacePot, *vq, *vr;
	double dq, dr, dis;
	char fileName[100], tmpChar[100];

	// dynamically allocate memory
	qSpacePot = (double *) calloc(NGRID, sizeof(double));
	rSpacePot = (double *) calloc(NGRID, sizeof(double));
	vq  = (double *) calloc(NGRID, sizeof(double));
	vr  = (double *) calloc(NGRID, sizeof(double));
	
	// grid point spacings
	dq = (QMAX-QMIN) / (double)(NGRID);
	dr = 0.02*TWOPI / ((double)(NGRID) * dq);
	  
	// q-space and real-space grid points 
	for (i = 0, dis = QMIN; i < NGRID; i++, dis += dq) vq[i] = dis;
	for (i = 0, dis = 0.0 ; i < NGRID; i++, dis += dr) vr[i] = dis;

	// Added new
	double volPerAtom = ( params.cellVolume / ((double)(params.nAtoms)) );
	printf("The cell volume = %.8f\n", params.cellVolume);
	printf("The volume per atom = %.8f\n", volPerAtom);
	fflush(stdout);

	// Calculate the q-space and real-space potentials along the grid points
	//int ina = 0;
	//for (k = 0; k < params.nBandStructures; k++) {
  	for (j = 0; j < params.nAtoms; j++) {
  		// calculate the q-space potential along the grid points
  		for (i = 0; i < NGRID; i++) {
  			qSpacePot[i] = calcPot(vq[i], atoms[j].ppParams);
  			qSpacePot[i] *= volPerAtom; 
  		}
  		// calculate the real-space potential along 1 to NGRID points
  		for (ir = 1; ir < NGRID; ir++) {
		    for (rSpacePot[ir] = 0.0, i = 0; i < NGRID; i++) {
		      rSpacePot[ir] += (vq[i] * sin(vq[i]*vr[ir]) * qSpacePot[i]);
		    } 
		    rSpacePot[ir] *= (dq * 4.0 * PIE / (pow(TWOPI, 3.0) * vr[ir]));
		}
		// calculate the real-space potential at grid point 0, sin(0)/(0) = 1
		ir = 0;
		for (rSpacePot[ir] = 0.0, i = 0; i < NGRID; i++) {
			rSpacePot[ir] += (vq[i] * vq[i] * qSpacePot[i]);
		}
		rSpacePot[ir] *= (dq * 4.0 * PIE / (pow(TWOPI, 3.0)));
		// print the real-space potential
		strcpy(fileName, "pot_");
		sprintf(tmpChar, "%s.dat", atoms[j].symbol);
		strcat(fileName, tmpChar);
		pf = fopen(fileName, "w"); 
		for (i = 0; i < NGRID; i++) {
			fprintf(pf, "%g %g\n", vr[i], rSpacePot[i]);
		}
		fclose(pf);
		// // increment the atom
  // 		ina++; 
  	}
  	//}

	// free dynamically allocated memory 
	free(qSpacePot); free(rSpacePot); 
	free(vq); free(vr);
	
	return;
}
  
/*****************************************************************************************/