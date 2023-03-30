/*****************************************************************************/
//
// This file deals with grid and gridPoint structures
//
/*****************************************************************************/

#include "qp.h"

/****************************************************************************/
// 

long initRSpaceGrid(grid *rSpaceGrid, gridPoint *rSpaceGP, lParams lPar) {
	FILE *pf;
	long iGrid, iX, iY, iZ, nAtoms;
	double x, y, z;

	// r-space grid point will be indexed by 0 (qSpace by 1)
	rSpaceGrid->index = 0;
	rSpaceGrid->gP = rSpaceGP;

	// Determine the dimensions of the r-space grid required 
	if (lPar.readInGrid) {
		readGridParFile(rSpaceGrid, "grid.par");
	}
	else if ( access("conf.par", F_OK) != -1 ) {
		pf = fopen("conf.par" , "r");
		fscanf(pf, "%ld", &nAtoms);
		determineRSpaceGridSize(rSpaceGrid, nAtoms, pf);
		fclose(pf);
	}
	else {
		printf("\n\nNo conf.par file detected in current working directory - the program is exiting!!!\n\n");
		fflush(stdout);
		exit(EXIT_FAILURE);
	}

	// Fill in the grid point structures for all grid points
	// The grid points go in the order of z->y->x (i.e. x changes the quickest)
	iGrid = 0;
	z = rSpaceGrid->minPos.z;  
	for (iZ = 0; iZ < rSpaceGrid->nGridPointsZ; iZ++) {
		y = rSpaceGrid->minPos.y;
		for (iY = 0; iY < rSpaceGrid->nGridPointsY; iY++) {
			x = rSpaceGrid->minPos.x;
			for (iX = 0; iX < rSpaceGrid->nGridPointsX; iX++) {
				rSpaceGP[iGrid].index = iGrid;
				rSpaceGP[iGrid].pos.z = z;
				rSpaceGP[iGrid].pos.y = y;
				rSpaceGP[iGrid].pos.x = x;
				rSpaceGP[iGrid].pos.mag = retVectorMagnitude(rSpaceGP[iGrid].pos);
				iGrid++;
				x += rSpaceGrid->stepSize.x;
			}
			y += rSpaceGrid->stepSize.y;
		}
		z += rSpaceGrid->stepSize.z;
	}

	return 0;
}

/****************************************************************************/
// 

long determineRSpaceGridSize(grid *rSpaceGrid, long nAtoms, FILE *confFilePointer) {
	long iAtom;
	double rx[nAtoms], ry[nAtoms], rz[nAtoms];
	atm_st atom[nAtoms];

	// Read in the atom positions to be used in determining the size of the r-space grid
	readConfFile(rx, ry, rz, atom, nAtoms, confFilePointer);

	// Box dimensions
	rSpaceGrid->maxPos.x = rint(0.5*retRangeOfDoubleArray(rx, nAtoms) + 5.0);
	rSpaceGrid->maxPos.y = rint(0.5*retRangeOfDoubleArray(ry, nAtoms) + 5.0);
	rSpaceGrid->maxPos.z = rint(0.5*retRangeOfDoubleArray(rz, nAtoms) + 5.0);
	rSpaceGrid->minPos.x = -rSpaceGrid->maxPos.x;
	rSpaceGrid->minPos.y = -rSpaceGrid->maxPos.y;
	rSpaceGrid->minPos.z = -rSpaceGrid->maxPos.z;
	rSpaceGrid->volume = 8.0*rSpaceGrid->maxPos.x*rSpaceGrid->maxPos.y*rSpaceGrid->maxPos.z;

	// dx, dy, dz, dr and dV=dx*dy*dz -> spacing between the grid points and the volume element
	rSpaceGrid->stepSize.x  = (rSpaceGrid->maxPos.x - rSpaceGrid->minPos.x) / (double)(rSpaceGrid->nGridPointsX);
	rSpaceGrid->stepSize.y  = (rSpaceGrid->maxPos.y - rSpaceGrid->minPos.y) / (double)(rSpaceGrid->nGridPointsY);
	rSpaceGrid->stepSize.z  = (rSpaceGrid->maxPos.z - rSpaceGrid->minPos.z) / (double)(rSpaceGrid->nGridPointsZ);
	rSpaceGrid->stepSize.mag = retVectorMagnitude(rSpaceGrid->stepSize);
	rSpaceGrid->dV = rSpaceGrid->stepSize.x*rSpaceGrid->stepSize.y*rSpaceGrid->stepSize.z;

	// Print grid related statistics
	fprintf(stdout, "Box (quadrant) dimensions: xMax = %.2f ymax = %.2f zMax = %.2f\n", 
						rSpaceGrid->maxPos.x, rSpaceGrid->maxPos.y, rSpaceGrid->maxPos.z);
	fprintf(stdout, "Box volume = %.2f\n", rSpaceGrid->volume); 
	fprintf(stdout, "Grid point spacing: dx = %.6f dy = %.6f dz = %.6f dr = %.6f\n", 
						rSpaceGrid->stepSize.x, rSpaceGrid->stepSize.y, rSpaceGrid->stepSize.z, rSpaceGrid->stepSize.mag);
 	fprintf(stdout, "Grid point volume, dV = %.6f\n", rSpaceGrid->dV);
 	fflush(stdout);

	return 0;
}

/****************************************************************************/
//

long retIndexOfClosestGridPoint(grid rSpaceGrid, vector targetPoint) {
	long iGrid, iClosestGP = -1;
	double smallestMag = 100000.;
	vector diff;

	for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
		diff = retSubtractedVectors(rSpaceGrid.gP[iGrid].pos, targetPoint);
		if (diff.mag < smallestMag) {
			iClosestGP = iGrid;
			smallestMag = diff.mag;
		}
	}

	return iClosestGP;
}

/****************************************************************************/
// 

long initQSpaceGrid(grid *qSpaceGrid, gridPoint *qSpaceGP, grid rSpaceGrid) {
	FILE *pf;
	long iGrid, iKx, iKy, iKz, iKyz;
	double *kx, *ky, *kz, OneOvernGridPoints;

	// q-space grid point will be indexed by 1 (rSpace by 0)
	qSpaceGrid->index = 0;
	qSpaceGrid->gP = qSpaceGP;
	qSpaceGrid->nGridPointsX = rSpaceGrid.nGridPointsX; 
	qSpaceGrid->nGridPointsY = rSpaceGrid.nGridPointsY;
	qSpaceGrid->nGridPointsZ = rSpaceGrid.nGridPointsZ;
	qSpaceGrid->nGridPoints = rSpaceGrid.nGridPoints;
	OneOvernGridPoints = 1.0 / ((double)qSpaceGrid->nGridPoints);

	// Determine the grid point spacing in q-space
	qSpaceGrid->stepSize.x = TWOPI / ((double)(rSpaceGrid.nGridPointsX) * rSpaceGrid.stepSize.x);
	qSpaceGrid->stepSize.y = TWOPI / ((double)(rSpaceGrid.nGridPointsY) * rSpaceGrid.stepSize.y);
	qSpaceGrid->stepSize.z = TWOPI / ((double)(rSpaceGrid.nGridPointsZ) * rSpaceGrid.stepSize.z);
	qSpaceGrid->stepSize.mag = retVectorMagnitude(qSpaceGrid->stepSize);
	qSpaceGrid->dV = qSpaceGrid->stepSize.x*qSpaceGrid->stepSize.y*qSpaceGrid->stepSize.z;

	// Allocate memory for temporary kx, ky and kz doubles
	if ((kx = (double *) calloc(qSpaceGrid->nGridPointsX, sizeof(double))) == NULL) memoryError("kx");
	if ((ky = (double *) calloc(qSpaceGrid->nGridPointsY, sizeof(double))) == NULL) memoryError("ky");
	if ((kz = (double *) calloc(qSpaceGrid->nGridPointsZ, sizeof(double))) == NULL) memoryError("kz");

	// Calculate kx, ky and kz (kx[0] = ky[0] = kz[0] = 0.0 as calloc was used)
	// OneOvernGridPoints not included here!
	for (iKx = 1; iKx <= (qSpaceGrid->nGridPointsX / 2); iKx++) {
		kx[iKx] = (double)(iKx) * qSpaceGrid->stepSize.x; // * OneOvernGridPoints;
		kx[qSpaceGrid->nGridPointsX - iKx] = -kx[iKx];
	}
	for (iKy = 1; iKy <= (qSpaceGrid->nGridPointsY / 2); iKy++) {
		ky[iKy] = (double)(iKy) * qSpaceGrid->stepSize.y; // * OneOvernGridPoints;
		ky[qSpaceGrid->nGridPointsY - iKy] = -ky[iKy];
	}
	for (iKz = 1; iKz <= (qSpaceGrid->nGridPointsZ / 2); iKz++) {
		kz[iKz] = (double)(iKz) * qSpaceGrid->stepSize.z; // * OneOvernGridPoints;
		kz[qSpaceGrid->nGridPointsZ - iKz] = -kz[iKz];
	}

	// Fill in the q-space grid point structures for all grid points
	// The grid points go in the order of z->y->x (i.e. x changes the quickest)
	for (iKz = 0; iKz < qSpaceGrid->nGridPointsZ; iKz++) {
		for (iKy = 0; iKy < qSpaceGrid->nGridPointsY; iKy++) {
			iKyz = qSpaceGrid->nGridPointsX * (qSpaceGrid->nGridPointsY * iKz + iKy);
			for (iKx = 0; iKx < qSpaceGrid->nGridPointsX; iKx++) {
				iGrid = iKyz + iKx;
				qSpaceGP[iGrid].index = iGrid;
				qSpaceGP[iGrid].pos.x = kx[iKx];
				qSpaceGP[iGrid].pos.y = ky[iKy];
				qSpaceGP[iGrid].pos.z = kz[iKz];
				qSpaceGP[iGrid].pos.mag = retVectorMagnitude(qSpaceGP[iGrid].pos);
			}
		}
	}

	// Fill in max and minimum position of the q-space grid
	qSpaceGrid->maxPos.x = -kx[(qSpaceGrid->nGridPointsX / 2)];
	qSpaceGrid->maxPos.y = -ky[(qSpaceGrid->nGridPointsY / 2)];
	qSpaceGrid->maxPos.z = -kz[(qSpaceGrid->nGridPointsZ / 2)];
	qSpaceGrid->minPos = retScaledVector(qSpaceGrid->maxPos, -1.0);
	qSpaceGrid->volume = qSpaceGrid->maxPos.x * qSpaceGrid->maxPos.y * qSpaceGrid->maxPos.z;

	// Print grid related statistics
	// TODO: turn this into a function as same code is used above in initRSpaceGrid
	fprintf(stdout, "Box (quadrant) dimensions: qxMax = %.2f qymax = %.2f qzMax = %.2f\n", 
						qSpaceGrid->maxPos.x, qSpaceGrid->maxPos.y, qSpaceGrid->maxPos.z);
	fprintf(stdout, "Q-space box volume = %.2f\n", qSpaceGrid->volume); 
	fprintf(stdout, "Grid point spacing: dqx = %.6f dqy = %.6f dqz = %.6f dqr = %.6f\n", 
						qSpaceGrid->stepSize.x, qSpaceGrid->stepSize.y, qSpaceGrid->stepSize.z, qSpaceGrid->stepSize.mag);
 	fprintf(stdout, "Grid point volume, dqV = %.6f\n", qSpaceGrid->dV);
 	fflush(stdout);

	// Free dynamically allocated memory
	free(kz); free(ky); free(kx);

	return 0;
}

/****************************************************************************/
