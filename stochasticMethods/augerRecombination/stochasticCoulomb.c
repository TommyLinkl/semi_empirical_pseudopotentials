/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "ar.h"

/*****************************************************************************/
// RrsZeta is ordered from iZeta then r then s

void calcRrsZetaMatrix(zomplex *RrsZeta, double *psiR, long nR, double *psiS, long nS, 
				zomplex *psiZeta, long nZeta, long nGridPoints, double dV) 
{
	long r, s, iZeta, iGrid;
	long nRS, iZetanGrid, iZetanRS, rnGrid, rnSiZetanRS, snGrid; // for optimization
	zomplex sum;

	nRS = nR*nS;
#pragma omp parallel for private(iZeta, r, s, iGrid, sum, iZetanGrid, iZetanRS, rnGrid, rnSiZetanRS, snGrid)
	for (iZeta = 0; iZeta < nZeta; iZeta++) {
		iZetanGrid = iZeta*nGridPoints;
		iZetanRS = iZeta*nRS;
		for (r = 0; r < nR; r++) {
			rnGrid = r*nGridPoints;
			rnSiZetanRS = r*nS + iZetanRS;
			for (s = 0; s < nS; s++) {
				snGrid = s*nGridPoints;
				sum.re = 0.0; sum.im = 0.0;
				for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
					sum.re += (psiS[snGrid + iGrid] * psiR[rnGrid + iGrid] * psiZeta[iZetanGrid + iGrid].re);
					sum.im += (psiS[snGrid + iGrid] * psiR[rnGrid + iGrid] * psiZeta[iZetanGrid + iGrid].im);
				}
				RrsZeta[rnSiZetanRS + s].re = sum.re*dV;
				RrsZeta[rnSiZetanRS + s].im = sum.im*dV;
			}
		}
	}

	return;
}

/*****************************************************************************/
// ChiaiZeta is ordered from iExc then Zeta then a then i
// flag = 0 for hole channel and flag = 1 for elec channel

void calcChiaiZetaMatrix(zomplex *ChiaiZeta, zomplex *RrsZeta, long nZeta, double *Cbs, lng_st ist, long elecHoleFlag) {
	long a, i, r, nA, nI, nR, iExc, iZeta; 
	long nAI, nZetaAI, anHoles, anR, anI, inR, iExcTHL, iExcnZetaAI, iZetanIR, iZetanAR, iZetanAI; // for optimization
	zomplex sum;

	// Useful constants to make code more readable
	const long nSO = nZeta;
	const long THL = ist.nNonIntExcitons;

	// Determine if electron or hole channel and set loops accordingly
	if (elecHoleFlag == 0) {      // RrsZeta = RijZeta
		nA = ist.nElecs; 
		nI = ist.itot;
		nR = ist.nHoles; // r = j and gets summed over
	}
	else if (elecHoleFlag == 1) { // RrsZeta = RabZeta
		nA = ist.atot; 
		nI = ist.nHoles;
		nR = ist.nElecs; // r = b and gets summed over
	}	
	else {
		printf("\n\nInvalid channel (0 (hole) or 1 (elec)) in function call to calcChiaiZetaMatrix!!!\n\n");
		exit(EXIT_FAILURE);
	}

	// Calculate ChiaiZeta for all permutations of: interacting excitonic states, stochastic orbitals, a and i states
	nAI = nA*nI;
	nZetaAI = nSO*nAI;
#pragma omp parallel for private(iExc, iZeta, a, i, r, sum, iExcTHL, iExcnZetaAI, iZetanAI, iZetanIR, iZetanAR, anHoles, anR, anI, inR)
	for (iExc = 0; iExc < ist.nIntExcitons; iExc++) {
		iExcTHL = iExc*THL;
		iExcnZetaAI = iExc*nZetaAI;
		for (iZeta = 0; iZeta < nSO; iZeta++) {
			iZetanAI = iZeta*nAI;
			iZetanIR = iZeta*nI*nR;
			iZetanAR = iZeta*nA*nR;
			for (a = 0; a < nA; a++) {
				anHoles = a*ist.nHoles;
				anR = a*nR;
				anI = a*nI;
				for (i = 0; i < nI; i++) {
					inR = i*nR;
					sum.re = 0.0; sum.im = 0.0;
					for (r = 0; r < nR; r++) {
						if (elecHoleFlag == 0) {      // hole channel sum over j of Caj*RijZeta
							sum.re += Cbs[iExcTHL + anR + r]*RrsZeta[iZetanIR + inR + r].re;
							sum.im += Cbs[iExcTHL + anR + r]*RrsZeta[iZetanIR + inR + r].im;
						}
						else if (elecHoleFlag == 1) { // elec channel sum over b of Cbi*RabZeta
							sum.re += Cbs[iExcTHL + r*ist.nHoles + i]*RrsZeta[iZetanAR + anR + r].re;
							sum.im += Cbs[iExcTHL + r*ist.nHoles + i]*RrsZeta[iZetanAR + anR + r].im;
						}
					}
					ChiaiZeta[iExcnZetaAI + iZetanAI + anI + i].re = sum.re;
					ChiaiZeta[iExcnZetaAI + iZetanAI + anI + i].im = sum.im;
				}
			}
		}
	}

	return;
}


/*****************************************************************************/
// TZeta is ordered from iExc then iZeta

void calcTZetaMatrix(zomplex *TZeta, zomplex *RckZeta, long nZeta, double *Cbs, lng_st ist) {
	long c, k, iExc, iZeta;
	long iExcnNonInt, iExcnZeta, iZetanNonInt, cnHoles; // for optimization 
	zomplex sum;

#pragma omp parallel for private(iExc, iZeta, c, k, sum, iExcnNonInt, iExcnZeta, iZetanNonInt, cnHoles)
	for (iExc = 0; iExc < ist.nIntExcitons; iExc++) {
		iExcnNonInt = iExc*ist.nNonIntExcitons;
		iExcnZeta = iExc*nZeta;
		for (iZeta = 0; iZeta < nZeta; iZeta++) {
			iZetanNonInt = iZeta*ist.nNonIntExcitons;
			sum.re = 0.0; sum.im = 0.0;
			for (c = 0; c < ist.nElecs; c++) {
				cnHoles = c*ist.nHoles;
				for (k = 0; k < ist.nHoles; k++) {
					sum.re += Cbs[iExcnNonInt + cnHoles + k]*RckZeta[iZetanNonInt + cnHoles + k].re;
					sum.im += Cbs[iExcnNonInt + cnHoles + k]*RckZeta[iZetanNonInt + cnHoles + k].im;
				}
			}
			TZeta[iExcnZeta + iZeta].re = sum.re; 
			TZeta[iExcnZeta + iZeta].im = sum.im;
		}
	}

	return;
}

/*****************************************************************************/
// Aai is ordered from iExc1 then iExc2 then a then i

void calcAaiMatrix(zomplex *Aai, zomplex *ChiaiZeta, zomplex *TZeta, long nZeta, long nIntExc, long nA, long nI) {
	long iZeta, iExc1, iExc2, a, i;
	long nAI, nZetaAI, iExc1nZetaAI, iExc2nZeta, iExc1nIntExciExc2nAI, anI, ianI;
	zomplex sum, tmp;

	nAI = nA*nI;
	nZetaAI = nZeta*nAI;
#pragma omp parallel for private(iExc1, iExc2, a, i, iZeta, sum, iExc1nZetaAI, iExc2nZeta, iExc1nIntExciExc2nAI, anI, ianI)
	for (iExc1 = 0; iExc1 < nIntExc; iExc1++) {
		iExc1nZetaAI = iExc1*nZetaAI;
		for (iExc2 = 0; iExc2 < nIntExc; iExc2++) {
			iExc2nZeta = iExc2*nZeta;
			iExc1nIntExciExc2nAI = ((iExc1*nIntExc + iExc2)*nAI);
			for (a = 0; a < nA; a++) {
				anI = a*nI;
				for (i = 0; i < nI; i++) {
					ianI = anI + i;
					sum.re = 0.0; sum.im = 0.0;
					for (iZeta = 0; iZeta < nZeta; iZeta++) {
						sum.re += ((ChiaiZeta[iExc1nZetaAI + iZeta*nAI + ianI].re * TZeta[iExc2nZeta + iZeta].re)
								  -(ChiaiZeta[iExc1nZetaAI + iZeta*nAI + ianI].im * TZeta[iExc2nZeta + iZeta].im));
						sum.im += ((ChiaiZeta[iExc1nZetaAI + iZeta*nAI + ianI].im * TZeta[iExc2nZeta + iZeta].re)
								  +(ChiaiZeta[iExc1nZetaAI + iZeta*nAI + ianI].re * TZeta[iExc2nZeta + iZeta].im));
					}
					Aai[iExc1nIntExciExc2nAI + anI + i].re += sum.re;
					Aai[iExc1nIntExciExc2nAI + anI + i].im += sum.im;
				}
			}
		}
	}

	return;
}

/*****************************************************************************/
// Calculate the stochastic orbitals (thetaZeta) used to approximate the Coulomb operator
// thetaZeta(r) = (1/TWOPI)^3*Integral(dk*sqrt(potq)*exp^(i*phi)*exp^(i*k.r) (backwards FT)
// phi is a random phase between 0 and TWOPI  
// potq is the FT of the Coulomb potential (scaled by 1/ist.nGirdPoints in init.c)
// FFT will be used to get the r-space orbitals thetaZeta(r)

long calcRandomCoulombStates(zomplex *thetaZeta, long nZeta, zomplex *potq, lng_st ist, par_st par, 
            fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi) 
{
  long iZeta, iGrid, iZetanGrid, idum = ist.seed;
  static long beenCalledBeforeFlag = 0;
  long nNegativePotQValues = 0; 
  double phi;
  double *sqrtPotk, sqrtdVolume_1 = 1.0/sqrt(par.dv);
  zomplex *tmpThetaZeta;

  // Allocate memory for Fourier transforms and sqrt of k-space Coulomb potential
  if ((tmpThetaZeta = (zomplex *) calloc(ist.nGridPoints, sizeof(zomplex))) == NULL) nerror("tmpThetaZeta");
  if ((sqrtPotk = (double *) calloc(ist.nGridPoints, sizeof(double))) == NULL) nerror("sqrtPotk");

  // Calculate the sqrt of the k-space Coulomb potential at each grid k-space grid point
  for (iGrid = 0; iGrid < ist.nGridPoints; iGrid++) {
  	if (potq[iGrid].re > 0.0) {
  		sqrtPotk[iGrid] = sqrt(potq[iGrid].re);
  	}
  	else if (potq[iGrid].re < (-1.0*EPS)) {
  		if (! beenCalledBeforeFlag) {
  			nNegativePotQValues++;
  		}
  		sqrtPotk[iGrid] = 0.0;
  	}
  }
  if (! beenCalledBeforeFlag && nNegativePotQValues) {
  	fprintf(stdout, "WARNING: potq has %ld negative values!\n", nNegativePotQValues); fflush(stdout);
  }

  // Create the random orbitals in k-space then FFT to get thetaZeta in r-space
  for (iZeta = 0; iZeta < nZeta; iZeta++) {
    // Generate a random thetaZeta in k-space
    for (iGrid = 0; iGrid < ist.nGridPoints; iGrid++) {
      phi = TWOPI*(ran_nrc(&idum)); // random phase between 0 and TWOPI
      tmpThetaZeta[iGrid].re = sqrtPotk[iGrid]*cos(phi);
      tmpThetaZeta[iGrid].im = sqrtPotk[iGrid]*sin(phi);
    }
    // FFT the k-space tmpThetaZeta to get the r-space thetaZeta
    memcpy(&fftwpsi[0], &tmpThetaZeta[0], ist.nGridPoints*sizeof(fftwpsi[0]));
    fftw_execute(planbw[0]);  
    memcpy(&tmpThetaZeta[0], &fftwpsi[0], ist.nGridPoints*sizeof(tmpThetaZeta[0])); 
	// Copy and scale to get final thetaZeta
    iZetanGrid = iZeta*ist.nGridPoints;
    for (iGrid = 0; iGrid < ist.nGridPoints; iGrid++) {
      thetaZeta[iZetanGrid + iGrid].re = tmpThetaZeta[iGrid].re*sqrtdVolume_1;
      thetaZeta[iZetanGrid + iGrid].im = tmpThetaZeta[iGrid].im*sqrtdVolume_1;
    }
  }

  // Set flag to 1 so that warning printing is only done the 1st time this function is called
  beenCalledBeforeFlag = 1;

  // Free dynamically allocated memory
  free(sqrtPotk);
  free(tmpThetaZeta);

  return idum;
}

/*****************************************************************************/
// A debugging/ testing function. Calculates the real space Coulomb potential
// from the set of stochastic Coulomb orbitals and calculates the error.
// Calculate <thetaZeta(r=0)*thetaZetaComplexConj(r)> = V(r) = 1/nSO * sum((tZ(r=0)*tZ*(r)))
// V(r) = (tZ(0).re + tZ(0).im) * (tZ*(r).re + tZ*(r).im)
// V(r) = (tZ(0).re + tZ(0).im) * (tZ(r).re - tZ(r).im)
// V(r).re = tZ(0).re*tZ(r).re + tZ(0).im*tZ(r).im
// V(r).im = tZ(0).im*tZ(r).re - tZ(0).re*tZ(r).im
// Then calculates the mean-squared error between V(r) and determinstic (FFT) of potq 
// Args: 
// 		thetaZeta = the ist.nStochOrbitals set of stochastic orbitals 
// 		potq = determinstic q-space potential that gets FT to create the deterministic r-space potential
// 		vx, vy, vz = real-space grid points 
// Output: 
// 		prints detRSpaceCoulombPot and stoRSpaceCoulombPot to rSpaceCoulombPot.dat

void calcStochRealSpaceCoulomb(zomplex *thetaZeta, zomplex *potq, double *vx, double *vy, double *vz, lng_st ist, par_st par, 
										fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwpsi) {
	FILE *pf;
	long i, ix, iy, iz, iyz, ixyz; 
	long iGrid, iOriginGridPoint, iZeta, iZetanGrid; 
	double x2, y2, z2, dr, boxLength;
	zomplex *detRSpaceCoulombPot, *stoRSpaceCoulombPot;
	vector rMinusRPrime;
	grid realSpaceGrid;

	// Write the beginning of the function
	fprintf(stdout, "\nStarting to calculate the error in the stochastic representation of the Coulomb potential\n");
	writeCurrentTime(stdout);
	fprintf(stdout, "The number of stochastic orbitals used = %ld\n", ist.nStochOrbitals);
	fflush(stdout);

	// Dynamically allocate memory
	if ((realSpaceGrid.gP = (gridPoint *) calloc(ist.nGridPoints, sizeof(gridPoint))) == NULL) nerror("gridPoint");
	if ((detRSpaceCoulombPot = (zomplex *) calloc(ist.nGridPoints, sizeof(zomplex))) == NULL) nerror("detRSpaceCoulombPot");
	if ((stoRSpaceCoulombPot = (zomplex *) calloc(ist.nGridPoints, sizeof(zomplex))) == NULL) nerror("stoRSpaceCoulombPot");

	// Fill in the index and position fields of the realSpaceGrid structure -> should be done
	// TODO: use realSpaceGrid throughout the program and pass it into this function
	realSpaceGrid.index = 1; 
	realSpaceGrid.stepSize.x = par.dx; realSpaceGrid.stepSize.y = par.dy; realSpaceGrid.stepSize.z = par.dz;
	realSpaceGrid.stepSize.mag = retVectorMagnitude(realSpaceGrid.stepSize);
	realSpaceGrid.dV = par.dx*par.dy*par.dz; // assumes rectangular grid with uniform spacing
	realSpaceGrid.nGridPoints = ist.nx*ist.ny*ist.nz;
	for (iz = 0; iz < ist.nz; iz++) {
		z2 = sqr(vz[iz]);
		for (iy = 0; iy < ist.ny; iy++) {
			y2 = sqr(vy[iy]);
			iyz = ist.nx * (ist.ny * iz + iy);
			for (ix = 0; ix < ist.nx; ix++) {
				x2 = sqr(vx[ix]);
				ixyz = iyz + ix;
				dr = sqrt(x2 + y2 + z2);
				realSpaceGrid.gP[ixyz].index = ixyz;
				realSpaceGrid.gP[ixyz].pos.x = vx[ix];
				realSpaceGrid.gP[ixyz].pos.y = vy[iy];
				realSpaceGrid.gP[ixyz].pos.z = vz[iz];
				realSpaceGrid.gP[ixyz].pos.mag = dr; 
			}
		}
	}

	// Calculate the deterministic real space potential
	// FFT the k-space Coulomb potential calculated in init.c to get the deterministic r-space Coulomb potential
    memcpy(&fftwpsi[0], &potq[0], ist.nGridPoints*sizeof(fftwpsi[0]));
    fftw_execute(planbw[0]);  
    memcpy(&detRSpaceCoulombPot[0], &fftwpsi[0], ist.nGridPoints*sizeof(detRSpaceCoulombPot[0]));

	// Determine index of the grid point corresponding to the origin, r = (0, 0, 0)
	for (iGrid = 0; iGrid < ist.nGridPoints; iGrid++) {
		if (realSpaceGrid.gP[iGrid].pos.mag < EPS) {
			iOriginGridPoint = iGrid;
			realSpaceGrid.iOrigin = iGrid;
			realSpaceGrid.origin = realSpaceGrid.gP[iGrid].pos;
			break;
		}
	}
	fprintf(stdout, "The index of the origin is = %ld\n", realSpaceGrid.iOrigin);
	fprintf(stdout, "Position of the origin is = % .4f % .4f % .4f\n", realSpaceGrid.gP[iOriginGridPoint].pos.x,  
							realSpaceGrid.gP[iOriginGridPoint].pos.y, realSpaceGrid.gP[iOriginGridPoint].pos.z);
	fflush(stdout);

	// Calculate <thetaZeta(r=0)*thetaZetaComplexConj(r)> = V(r) = 1/nSO * sum((tZ(r=0)*tZ*(r)))
	// V(r) = (tZ(r=0).re + tZ(r=0).im) * (tZ*(r).re + tZ*(r).im)
	for (iZeta = 0; iZeta < ist.nStochOrbitals; iZeta++) {
		iZetanGrid = iZeta*ist.nGridPoints; // for optimization
		for (iGrid = 0; iGrid < ist.nGridPoints; iGrid++) {
			// V(r).re = tZ(0).re*tZ(r).re + tZ(0).im*tZ(r).im
			stoRSpaceCoulombPot[iGrid].re += thetaZeta[iZetanGrid + iOriginGridPoint].re*thetaZeta[iZetanGrid + iGrid].re +
											thetaZeta[iZetanGrid + iOriginGridPoint].im*thetaZeta[iZetanGrid + iGrid].im;
			// V(r).im = tZ(0).im*tZ(r).re - tZ(0).re*tZ(r).im
			stoRSpaceCoulombPot[iGrid].im += thetaZeta[iZetanGrid + iOriginGridPoint].im*thetaZeta[iZetanGrid + iGrid].re -
											thetaZeta[iZetanGrid + iOriginGridPoint].re*thetaZeta[iZetanGrid + iGrid].im;							
		}
	}

	// Normalize the stochastic representation of the Coulomb potential by the number of stochastic orbitals
	for (iGrid = 0; iGrid < ist.nGridPoints; iGrid++) {
		stoRSpaceCoulombPot[iGrid].re /= ((double)(ist.nStochOrbitals));
		stoRSpaceCoulombPot[iGrid].im /= ((double)(ist.nStochOrbitals));
	}

	// Write a file that contains the deterministic and stochastic representations of the Coulomb potential
	pf = fopen("rSpaceCoulombPot.dat", "w");
	for (iGrid = 0; iGrid < ist.nGridPoints; iGrid++) {
		fprintf(pf, "%ld % .4f % .4f % .4f %.4f ", realSpaceGrid.gP[iGrid].index, realSpaceGrid.gP[iGrid].pos.x, 
					realSpaceGrid.gP[iGrid].pos.y, realSpaceGrid.gP[iGrid].pos.z, realSpaceGrid.gP[iGrid].pos.mag);
		fprintf(pf, "% .12f % .12f ", detRSpaceCoulombPot[iGrid].re, stoRSpaceCoulombPot[iGrid].re); 
		fprintf(pf, "% .12f % .12f\n", detRSpaceCoulombPot[iGrid].im, stoRSpaceCoulombPot[iGrid].im);		
	}
	fclose(pf);

	// Print function ending
	fprintf(stdout, "Finished calculating the stochastic representation of the Coulomb potential\n");
	writeCurrentTime(stdout);
	fflush(stdout);

	// Free dynamically allocated memory
	free(detRSpaceCoulombPot);
	free(stoRSpaceCoulombPot);
	free(realSpaceGrid.gP);

	return; 
}

/*****************************************************************************/
//

void printCoulombAndRMatricesConvergenceInfo(zomplex *RckZeta, zomplex *RabZeta, lng_st ist) {
	FILE *pf;
	long iZeta, c, k, a, b, zckIndex, zabIndex;
	double magSum;
	zomplex sum;

	// Check convergence of 1/nS0 * sum (<tZ|ck>)
	pf = fopen("RckConvergence.dat", "w");		
	for (c = 0; c < ist.nElecs; c++) {
		for (k = 0; k < ist.nHoles; k++) {
			sum.re = 0.0; sum.im = 0.0; magSum = 0.0;
			for (iZeta = 0; iZeta < ist.nStochOrbitals; iZeta++) {
				zckIndex = iZeta*ist.nHoles*ist.nElecs + c*ist.nHoles + k;
				sum.re += RckZeta[zckIndex].re;
				sum.im += RckZeta[zckIndex].im;
				magSum += sqr(RckZeta[zckIndex].re)+sqr(RckZeta[zckIndex].im);
				if ((c%5)==0 && (k%5)==0 && (iZeta%5)==0) {
					fprintf(pf, "%ld %ld %ld % g % g % g % g\n", c, k, iZeta+1, 
						sum.re/(double)(iZeta+1), sum.im/(double)(iZeta+1), 
						(sum.re*sum.re+sum.im*sum.im)/(double)(iZeta+1)/(double)(iZeta+1),
						magSum/(double)(iZeta+1));
				}
			}
		}
	}
	fclose(pf);

	// Check convergence of 1/nS0 * sum (<ab|tZ>)
	pf = fopen("RabConvergence.dat", "w");		
	for (a = 0; a < ist.atot; a++) {
		for (b = 0; b < ist.nElecs; b++) {
			sum.re = 0.0; sum.im = 0.0; magSum = 0.0;
			for (iZeta = 0; iZeta < ist.nStochOrbitals; iZeta++) {
				zabIndex = iZeta*ist.atot*ist.nElecs + a*ist.nElecs + b;
				sum.re += RabZeta[zabIndex].re;
				sum.im += RabZeta[zabIndex].im;
				magSum += sqr(RabZeta[zabIndex].re)+sqr(RabZeta[zabIndex].im);
				if ((a%5)==0 && (b%5)==0 && (iZeta%5)==0) {
					fprintf(pf, "%ld %ld %ld % g % g % g % g\n", a, b, iZeta+1, 
						sum.re/(double)(iZeta+1), sum.im/(double)(iZeta+1), 
						(sum.re*sum.re+sum.im*sum.im)/(double)(iZeta+1)/(double)(iZeta+1),
						magSum/(double)(iZeta+1));
				}
			}
		}
	}
	fclose(pf);

	// Check convergence of <ab|V|ck> = 1/nS0 * sum (<ab|tZ> <tZ|ck>)
	pf = fopen("vabck_stochastic.dat", "w");		
	for (a = 0; a < ist.atot; a++) {
		for (b = 0; b < ist.nElecs; b++) {
    		for (c = 0; c < ist.nElecs; c++) {
				for (k = 0; k < ist.nHoles; k++) {
					sum.re = 0.0; sum.im = 0.0;
					for (iZeta = 0; iZeta < ist.nStochOrbitals; iZeta++) {
						zckIndex = iZeta*ist.nHoles*ist.nElecs + c*ist.nHoles + k;
						zabIndex = iZeta*ist.atot*ist.nElecs + a*ist.nElecs + b;
						sum.re += RabZeta[zabIndex].re*RckZeta[zckIndex].re 
								- RabZeta[zabIndex].im*RckZeta[zckIndex].im;
						sum.im += RabZeta[zabIndex].re*RckZeta[zckIndex].im
								+ RabZeta[zabIndex].im*RckZeta[zckIndex].re;
					}
					sum.re /= (double)(ist.nStochOrbitals);
					sum.im /= (double)(ist.nStochOrbitals);
					fprintf(pf, "%ld %ld %ld %ld % g % g % g\n", a, b, c, k, 
								sum.re, sum.im, sum.re*sum.re+sum.im*sum.im);
    			}
			}
		}
	}
	fclose(pf);

	return;
}

/*****************************************************************************/
