/*****************************************************************************/
//
// This file calculates absorption related properties for excitonic states
//
/*****************************************************************************/

#include "qp.h"

/*****************************************************************************/
//

long calcAbsorptionProperties(intTwoQPState *intTwoQP, nonintTwoQPState *nonintTwoQP, grid rSpaceGrid, grid qSpaceGrid,
					dParams dPar, lParams lPar, fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwPsi) {
	FILE *pf;
	long i, iIntTwoQP, iNonintTwoQP;
	double *rotStrength, *intRotStrength;
	vector *elecDipoleMoment, *magDipoleMoment;
	vector *intElecDipMoment, *intMagDipMoment;

	// Write beginning of function
	writeSeparation(stdout);
	writeCurrentTime(stdout);
	fprintf(stdout, "Beginnning the calculation of the absorption properties\n\n");
	fflush(stdout);

	// Dynamically allocate memory
	if ((rotStrength = (double *) calloc(lPar.nNonintTwoQPStates, sizeof(double))) == NULL) memoryError("rotStrength");
	if ((magDipoleMoment = (vector *) calloc(lPar.nNonintTwoQPStates, sizeof(vector))) == NULL) memoryError("magDipoleMoment");
	if ((elecDipoleMoment = (vector *) calloc(lPar.nNonintTwoQPStates, sizeof(vector))) == NULL) memoryError("elecDipoleMoment");
	if ((intRotStrength = (double *) calloc(lPar.nIntTwoQPStates, sizeof(double))) == NULL) memoryError("intRotStrength");
	if ((intMagDipMoment = (vector *) calloc(lPar.nIntTwoQPStates, sizeof(vector))) == NULL) memoryError("intMagDipMoment");
	if ((intElecDipMoment = (vector *) calloc(lPar.nIntTwoQPStates, sizeof(vector))) == NULL) memoryError("intElecDipMoment");

	// Calculate absorption spectrum under the electric dipole approximation
	calcNonintElecDipoleMoments(elecDipoleMoment, nonintTwoQP, rSpaceGrid, lPar.nNonintTwoQPStates);
	calcIntElecDipoleMoments(intElecDipMoment, elecDipoleMoment, intTwoQP, lPar.nNonintTwoQPStates, lPar.nIntTwoQPStates);
	calcRadiativeLifetimesElecDipoleApprox(intElecDipMoment, elecDipoleMoment, intTwoQP, nonintTwoQP, dPar, lPar);

	// Calculate the magnetic dipole based absorption properties 
	// TODO: write below functions
	calcNonintMagneticMoments(magDipoleMoment, nonintTwoQP, rSpaceGrid, qSpaceGrid, lPar.nNonintTwoQPStates, planfw, planbw, fftwPsi);
	calcIntMagneticMoments(intMagDipMoment, magDipoleMoment, intTwoQP, lPar.nNonintTwoQPStates, lPar.nIntTwoQPStates);
	// calcRadiativeLifetimesMagDipole();

	// Calculate the rotational strength
	// calcNonintRotationalStrengths();
	// calcIntRotationalStrengths();

	// TODO: move into separate function or calcIntElecDipoleMoments at least
	writeIntTwoQPVectorProperty(intElecDipMoment, intTwoQP, lPar.nIntTwoQPStates, "intElecDipoleMoments.dat");

	if (! strcmp(lPar.calcType, "sp-elec-hole")) {
		for (i = 0; i < lPar.nNonintTwoQPStates; i++) {
			nonintTwoQP[i].spinZ = retSpinZProjNonintTwoQPState(nonintTwoQP[i]);
			nonintTwoQP[i].spin = retTotalSpinNonintTwoQPState(nonintTwoQP[i]);
			intTwoQP[i].spinZ = retSpinZProjIntTwoQPState(intTwoQP[i], lPar.nNonintTwoQPStates);
			intTwoQP[i].spin = retTotalSpinIntTwoQPState(intTwoQP[i], lPar.nNonintTwoQPStates);
		}
		writeSpinPolarizedNonintTwoQPVectorProperty(elecDipoleMoment, nonintTwoQP, lPar.nNonintTwoQPStates, "spNonintElecDipoleMoments.dat");
		writeSpinPolarizedIntTwoQPVectorProperty(intElecDipMoment, intTwoQP, lPar.nIntTwoQPStates, "spIntElecDipoleMoments.dat");
	}

	// Free dynamically allocated memory
	free(rotStrength); free(intRotStrength);
	free(magDipoleMoment); free(intMagDipMoment);
	free(elecDipoleMoment); free(intElecDipMoment);

	// Write ending of function
	fprintf(stdout, "\nFinished the calculation of the absorption properties\n");
	writeCurrentTime(stdout);
	fflush(stdout);

	return 0;
}

/*****************************************************************************/
//

long calcIntMagneticMoments(vector *intMagDipMoment, vector *magDipoleMoment, intTwoQPState *intTwoQP, 
								long nNonintStates, long nIntStates) {
	FILE *pf;
	long iIntState, iNonIntState;
	vector tmpVector, sum;

//#pragma omp parallel for private(iIntState, sum, iNonIntState, tmpVector)
	for (iIntState = 0; iIntState < nIntStates; iIntState++) {
		sum = retZeroVector();
		for (iNonIntState = 0; iNonIntState < nNonintStates; iNonIntState++) {
			tmpVector = retScaledVector(magDipoleMoment[iNonIntState], intTwoQP[iIntState].Crs[iNonIntState]);
			sum = retAddedVectors(sum, tmpVector);
		}
		intMagDipMoment[iIntState] = sum;
	}

	return 0;
}

/*****************************************************************************/
//

long calcNonintAbsorptionProperties(nonintTwoQPState *nonintTwoQP, grid rSpaceGrid, grid qSpaceGrid, dParams dPar, lParams lPar,
									fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwPsi) {
	FILE *pf;
	long iNonintTwoQP;
	double *rotStrength;
	vector *elecDipoleMoment, *magDipoleMoment;

	// Write beginning of function
	writeSeparation(stdout); writeCurrentTime(stdout);
	fprintf(stdout, "Beginnning the calculation of the noninteracting absorption properties\n\n");	fflush(stdout);

	// Dynamically allocate memory
	if ((rotStrength = (double *) calloc(lPar.nNonintTwoQPStates, sizeof(double))) == NULL) memoryError("rotStrength");
	if ((magDipoleMoment = (vector *) calloc(lPar.nNonintTwoQPStates, sizeof(vector))) == NULL) memoryError("magDipoleMoment");
	if ((elecDipoleMoment = (vector *) calloc(lPar.nNonintTwoQPStates, sizeof(vector))) == NULL) memoryError("elecDipoleMoment");

	// Calculate absorption spectrum under the electric dipole approximation
	calcNonintElecDipoleMoments(elecDipoleMoment, nonintTwoQP, rSpaceGrid, lPar.nNonintTwoQPStates);
	// TODO: write below function
	//calcNonintRadiativeLifetimesElecDipoleApprox(elecDipoleMoment, nonintTwoQP, dPar, lPar);

	// Calculate the magnetic dipole based absorption properties 
	// TODO: write / finish below functions
	calcNonintMagneticMoments(magDipoleMoment, nonintTwoQP, rSpaceGrid, qSpaceGrid, lPar.nNonintTwoQPStates, planfw, planbw, fftwPsi);
	// calcNonintRadiativeLifetimesMagDipole();

	// Calculate the rotational strength
	// TODO: write below function
	// calcNonintRotationalStrengths();

	// Print out a file with the noninteracting electric dipole matrix elements
	writeNonintTwoQPVectorProperty(elecDipoleMoment, nonintTwoQP, lPar.nNonintTwoQPStates, "nonintElecDipoleMoments.dat");
	
	// Free dynamically allocated memory
	free(rotStrength); free(magDipoleMoment); free(elecDipoleMoment); 

	// Write ending of function
	fprintf(stdout, "\nFinished the calculation of the noninteracting absorption properties\n");
	writeCurrentTime(stdout);	fflush(stdout);

	return 0;
}

/*****************************************************************************/
//

long calcNonintMagneticMoments(vector *magDipoleMoment, nonintTwoQPState *nonintTwoQP, grid rSpaceGrid, grid qSpaceGrid, 
								long nStates, fftw_plan_loc *planfw, fftw_plan_loc *planbw, fftw_complex *fftwPsi) {
	FILE *pf;
	long iState, iR, iS, iGrid;
	double tmp, *psiR, *psiS, *psiSdx, *psiSdy, *psiSdz;
	zomplex *cPsiS;
	vector tmpVector, sum, *firstDerivativeOfPsiS;

	// Useful constants
	const long nGridPoints = rSpaceGrid.nGridPoints;
	const double dV = rSpaceGrid.dV;

	// Dynamically allocate memory
	if ((psiSdx = (double *) calloc(nGridPoints, sizeof(double))) == NULL) memoryError("psiSdx");
	if ((psiSdy = (double *) calloc(nGridPoints, sizeof(double))) == NULL) memoryError("psiSdy");
	if ((psiSdz = (double *) calloc(nGridPoints, sizeof(double))) == NULL) memoryError("psiSdz");
	if ((cPsiS = (zomplex *) calloc(nGridPoints, sizeof(zomplex))) == NULL) memoryError("cPsiS");
	if ((firstDerivativeOfPsiS = (vector *) calloc(nGridPoints, sizeof(vector))) == NULL) memoryError("firstDerivativeOfPsiS");

	// Calculate <>
	for (iState = 0; iState < nStates; iState++) {
		psiR = nonintTwoQP[iState].qp2->psi;
		psiS = nonintTwoQP[iState].qp1->psi;
		sum = retZeroVector();
		// FFT psiS and store first derivative in firstDerivativeOfPsiS = (dPsiS/dx, dPsiS/dy, dPsiS/dz)
		// TODO: 
		//
		//for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
		//	tmp = psiR[iGrid]*psiS[iGrid];
		//	tmpVector = retScaledVector(rSpaceGrid.gP[iGrid].pos, tmp);
		//	sum = retAddedVectors(sum, tmpVector);
		//}
		//sum = retScaledVector(sum, dV);
		magDipoleMoment[iState] = sum;
	}

	// Free dynamically allocated memory
	free(firstDerivativeOfPsiS); free(cPsiS);
	free(psiSdz); free(psiSdy); free(psiSdx);

	return 0;
}

/*****************************************************************************/
//

long calcRadiativeLifetimesElecDipoleApprox(vector *intElecDipMoment, vector *elecDipoleMoment, intTwoQPState *intTwoQP,
										nonintTwoQPState *nonintTwoQP, dParams dPar, lParams lPar) {
	FILE *pf;
	long iIntTwoQP, iNonintTwoQP, i, iBin, index, nGaussianBins, nEmissionBins = 2000;
	double partitionFunction, lowestEnergy;
	double tmp, beta, temperature, totalRadiativeLifetime, totalAbsorptionCrossSection, boltzWeightedRadiativeRate, boltzWeightedRadiativeLifetime;
	double energy, scaling, dE, speedOfLight = 137.0;
	double emissionEnergy[nEmissionBins], emissionIntensity[nEmissionBins], absorptionIntensity[nEmissionBins];
	double emissionBinWidth = 0.0025/AUTOEV; // 5 meV bin width
	double minEmissionEnergy = 0.5/AUTOEV; // 500 meV minimum emission energy
	//double maxEmissionEnergy = (minEmissionEnergy + emissionBinWidth*((double)(nEmissionBins)));
	double emissionBroadeningSigma = dPar.energyLevelSigma; // sigma in gaussian for emission spectra broadening
	double emissionGuassianPrefactor = ( 1.0 / (emissionBroadeningSigma*sqrt(2.0*PIE)) ); 
	vector radLifetime; 

	// Set constants
	temperature = dPar.electronicTemperature;
	scaling = 3.0*cube(speedOfLight)*AUTONS/4.0;
	beta = AUTOEV/(KB*temperature);
	nGaussianBins = 2*(long)( 4.0*2.0*emissionBroadeningSigma / emissionBinWidth);

	// Fill in emission energy bins
	emissionEnergy[0] = minEmissionEnergy;
	for (iBin = 1; iBin < nEmissionBins; iBin++) {
		emissionEnergy[iBin] = ( emissionEnergy[iBin-1] + emissionBinWidth );
	}
	zeroDoubleArray(emissionIntensity, nEmissionBins);
	zeroDoubleArray(absorptionIntensity, nEmissionBins);

	// Noninteracting elec-hole pairs
	pf = fopen("nonintRadiativeLifetimes.dat", "w");
	lowestEnergy = nonintTwoQP[0].energy;
	partitionFunction = 0.0;
	boltzWeightedRadiativeRate = 0.0;
	for (iNonintTwoQP = 0; iNonintTwoQP < lPar.nNonintTwoQPStates; iNonintTwoQP++) {
		// Radiative lifetime = 3c^3/(4*omega^3*|u_12|^2) = scaling/(energy^3*elecDipoleMoment.mag^2)
		energy = nonintTwoQP[iNonintTwoQP].energy;
		tmp = scaling / cube(energy);
		radLifetime.x = tmp / sqr(elecDipoleMoment[iNonintTwoQP].x);
		radLifetime.y = tmp / sqr(elecDipoleMoment[iNonintTwoQP].y);
		radLifetime.z = tmp / sqr(elecDipoleMoment[iNonintTwoQP].z);
		totalRadiativeLifetime = tmp / sqr(elecDipoleMoment[iNonintTwoQP].mag);
                // added by Dipti: emission dependsd on energy^3, absorption on energy
                totalAbsorptionCrossSection = tmp * sqr(energy) / sqr(elecDipoleMoment[iNonintTwoQP].mag);
		fprintf(pf, "%5ld %3ld %3ld % .8f % .8f  %.2f  %4.5f  %g  %g  %g\n", iNonintTwoQP, nonintTwoQP[iNonintTwoQP].qp1->index,
						nonintTwoQP[iNonintTwoQP].qp2->index, energy, energy*AUTOEV, AUTOLAMBDANM/energy, 
						totalRadiativeLifetime, radLifetime.x, radLifetime.y, radLifetime.z);
		tmp = exp(-beta * (energy - lowestEnergy));
		partitionFunction += tmp;
		boltzWeightedRadiativeRate += tmp/totalRadiativeLifetime;
		// Add boltzmann weighted emission intensity/rate to correct emissionIntensity bin, Gaussian broadened
		iBin = (long)( (energy - minEmissionEnergy) / emissionBinWidth ); // this is the center bin
		for (i = (iBin - nGaussianBins/2); i <= (iBin + nGaussianBins/2); i++) {
			if (i < nEmissionBins) {
				dE = ( emissionEnergy[iBin] - emissionEnergy[i] );
				absorptionIntensity[i] += (1.0/totalAbsorptionCrossSection * emissionGuassianPrefactor * exp(-0.5*sqr(dE/emissionBroadeningSigma))); 
				emissionIntensity[i]   += (tmp/totalRadiativeLifetime * emissionGuassianPrefactor * exp(-0.5*sqr(dE/emissionBroadeningSigma)));
			}
		}
	}
	boltzWeightedRadiativeRate /= partitionFunction;
	boltzWeightedRadiativeLifetime = 1.0/boltzWeightedRadiativeRate;
	fclose(pf);	

	// Print Boltzmann weighted statistics
	fprintf(stdout, "Noninteracting elec-hole pairs Boltzmann weighted radiative recombination statistics:\n");
	fprintf(stdout, "Temperature = %.1f\n", temperature);
	fprintf(stdout, "Partition function = %.3f\n", partitionFunction);
	fprintf(stdout, "Electric dipole radiative rate = %.6f ns^-1\n", boltzWeightedRadiativeRate); 
	fprintf(stdout, "Electric dipole radiative lifetime = %.6f ns\n", boltzWeightedRadiativeLifetime); fflush(stdout);

	// Normalize the emission spectra and print the absoprtion and emission spectra
	normalizeDoubleArray(emissionIntensity, nEmissionBins);
	pf = fopen("nonintAbsorptionEmissionSpectra.dat", "w");
	fprintf(pf, "#i        eV        nm        absoprtion           emission\n");
	for (iBin = 0; iBin < nEmissionBins; iBin++) {
		fprintf(pf, "%4ld  %.8f %8.2f %.16f  %.16f\n", iBin, emissionEnergy[iBin]*AUTOEV, AUTOLAMBDANM/emissionEnergy[iBin],
														absorptionIntensity[iBin], emissionIntensity[iBin]);
	}
	fclose(pf);

	// Interacting elec-hole pairs
	zeroDoubleArray(emissionIntensity, nEmissionBins);
	zeroDoubleArray(absorptionIntensity, nEmissionBins);
	pf = fopen("intRadiativeLifetimes.dat", "w");
	fprintf(pf, "#iwoQD  eTwoQP[au]   eTwoQP[eV]   eTwoQP[nm]    totalRadiativeLifetime    radLifetime.x   radLifetime.y    radLifetime.z\n");
	lowestEnergy = intTwoQP[0].energy;
	partitionFunction = 0.0;
	boltzWeightedRadiativeRate = 0.0;
        printf("hi %ld\n", lPar.nIntTwoQPStates); 
	for (iIntTwoQP = 0; iIntTwoQP < lPar.nIntTwoQPStates; iIntTwoQP++) {
		// Radiative lifetime = 3c^3/(4*omega^3*|u_12|^2) = scaling/(energy^3*intElecDipMoment.mag^2)
		energy = intTwoQP[iIntTwoQP].energy;
		tmp = scaling / cube(energy);
		radLifetime.x = tmp / sqr(intElecDipMoment[iIntTwoQP].x);
		radLifetime.y = tmp / sqr(intElecDipMoment[iIntTwoQP].y);
		radLifetime.z = tmp / sqr(intElecDipMoment[iIntTwoQP].z);
		totalRadiativeLifetime = tmp / sqr(intElecDipMoment[iIntTwoQP].mag);
		totalAbsorptionCrossSection = tmp * sqr(energy) / sqr(intElecDipMoment[iIntTwoQP].mag);
		fprintf(pf, "%5ld % .8f % .8f %.2f  %4.5f  %g  %g  %g\n", iIntTwoQP, energy, energy*AUTOEV, AUTOLAMBDANM/energy, totalRadiativeLifetime, radLifetime.x, radLifetime.y, radLifetime.z);
		tmp = exp(-beta * (energy - lowestEnergy));
		partitionFunction += tmp;
		boltzWeightedRadiativeRate += tmp/totalRadiativeLifetime;
		// Add boltzmann weighted emission intensity/rate to correct emissionIntensity bin, Gaussian broadened
		iBin = (long)( (energy - minEmissionEnergy) / emissionBinWidth ); // this is the center bin
		for (i = (iBin - nGaussianBins/2); i <= (iBin + nGaussianBins/2); i++) {
			if (i < nEmissionBins) {
				dE = ( emissionEnergy[iBin] - emissionEnergy[i] );
				// TODO: check accuracy of cube(energy) from absorption intensity calculation
				absorptionIntensity[i] += (1.0/totalAbsorptionCrossSection * emissionGuassianPrefactor * exp(-0.5*sqr(dE/emissionBroadeningSigma))); 
				// emissionIntensity[i]   += (tmp/totalRadiativeLifetime * emissionGuassianPrefactor * exp(-0.5*sqr(dE/emissionBroadeningSigma)));
				emissionIntensity[i]   += (1.0/totalRadiativeLifetime * emissionGuassianPrefactor * exp(-0.5*sqr(dE/emissionBroadeningSigma)));
			}
		}
	}
	boltzWeightedRadiativeRate /= partitionFunction;
	boltzWeightedRadiativeLifetime = 1.0/boltzWeightedRadiativeRate;
	fclose(pf);	

	// Print Boltzmann weighted statistics
	fprintf(stdout, "\nInteracting elec-hole pairs Boltzmann weighted radiative recombination statistics:\n");
	fprintf(stdout, "Temperature = %.1f\n", temperature);
	fprintf(stdout, "Partition function = %.3f\n", partitionFunction);
	fprintf(stdout, "Electric dipole radiative rate = %.6f ns^-1\n", boltzWeightedRadiativeRate); 
	fprintf(stdout, "Electric dipole radiative lifetime = %.6f ns\n\n", boltzWeightedRadiativeLifetime); fflush(stdout);
	
	// Normalize the emission spectra and print the absoprtion and emission spectra
	normalizeDoubleArray(emissionIntensity, nEmissionBins);
	pf = fopen("intAbsorptionEmissionSpectra.dat", "w");
	fprintf(pf, "#i        eV        nm        absoprtion           emission\n");
	for (iBin = 0; iBin < nEmissionBins; iBin++) {
		fprintf(pf, "%4ld  %.8f %8.2f %.16f  %.16f\n", iBin, emissionEnergy[iBin]*AUTOEV, AUTOLAMBDANM/emissionEnergy[iBin],
														absorptionIntensity[iBin], emissionIntensity[iBin]);
	}
	fclose(pf);

	// Write temperature and broadening used in absorption / emission spectra calculations
	fprintf(stdout, "Electronic temperature in emission spectrum = %.1f K\n", temperature);
	fprintf(stdout, "Sigma used in Gaussian broadening of absorption and emission spectra = %.1f meV\n", emissionBroadeningSigma*AUTOMEV);
	fflush(stdout);

	return 0;
}

/*****************************************************************************/
//

long calcIntElecDipoleMoments(vector *intElecDipMoment, vector *elecDipoleMoment, intTwoQPState *intTwoQP, 
								long nNonintStates, long nIntStates) {
	FILE *pf;
	long iIntState, iNonIntState;
	vector tmpVector, sum;

//#pragma omp parallel for private(iIntState, sum, iNonIntState, tmpVector)
	for (iIntState = 0; iIntState < nIntStates; iIntState++) {
		sum = retZeroVector();
		for (iNonIntState = 0; iNonIntState < nNonintStates; iNonIntState++) {
			tmpVector = retScaledVector(elecDipoleMoment[iNonIntState], intTwoQP[iIntState].Crs[iNonIntState]);
			sum = retAddedVectors(sum, tmpVector);
		}
		intElecDipMoment[iIntState] = sum;
	}

	return 0;
}

/*****************************************************************************/
//

long calcNonintElecDipoleMoments(vector *elecDipoleMoment, nonintTwoQPState *nonintTwoQP, grid rSpaceGrid, long nStates) {
	FILE *pf;
	long iState, iGrid;
	double tmp, *psiR, *psiS;
	vector tmpVector, sum;

	// Useful constants
	const long nGridPoints = rSpaceGrid.nGridPoints;
	const double dV = rSpaceGrid.dV;

//#pragma omp parallel for private(iState, psiR, psiS, sum, iGrid, tmp, tmpVector)
	for (iState = 0; iState < nStates; iState++) {
		psiR = nonintTwoQP[iState].qp2->psi;
		psiS = nonintTwoQP[iState].qp1->psi;
		sum = retZeroVector();
		for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
			tmp = psiR[iGrid]*psiS[iGrid];
			tmpVector = retScaledVector(rSpaceGrid.gP[iGrid].pos, tmp);
			sum = retAddedVectors(sum, tmpVector);
		}
		sum = retScaledVector(sum, dV);
		elecDipoleMoment[iState] = sum;
	}

	return 0;
}

/*****************************************************************************/




