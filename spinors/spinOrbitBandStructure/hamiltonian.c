/****************************************************************************/
/* this file calculates the all the matrix elements for the Hamiltonian matrix */

#include "fit.h"
#define calcBessel(x,eps)   ((x) < (eps) ? 0.0 : (sin((x)) / ((x)*(x)) - cos((x)) / (x))) 

/****************************************************************************/

// simply adds up the kinetic energy, local potential, and spin-orbit to get the hamiltonian matrix
void calcHamiltonianMatrix(complexnumber *hamiltonian, complexnumber *kineticEnergyMat, complexnumber *localPotential, complexnumber *spinOrbitMatrix, param params) {
  int i, j;
  int n = params.nBasisVectors;

  zeroComplexMatrix(hamiltonian, 2*n, 2*n);

  // add kinetic and local potential energies to the Hamiltonian (block-diagonal)
  // delta function on the spin results terms added in a block-diagonal structure 
  for (i = 0; i < n; i++) { 
    for (j = 0; j < n; j++) {
      // the kinetic energy matrix is completely real
      // up-up (top left) block of the matrix
      hamiltonian[2*i*n+j].re += localPotential[2*i*n+j].re + kineticEnergyMat[2*n*i+j].re;
      hamiltonian[2*i*n+j].im += localPotential[2*i*n+j].im;
      // down-down (bottom right) block of the matrix
      hamiltonian[2*n*n+n+2*n*i+j].re += localPotential[2*n*n+n+2*n*i+j].re + kineticEnergyMat[2*n*n+n+2*n*i+j].re;
      hamiltonian[2*n*n+n+2*n*i+j].im += localPotential[2*n*n+n+2*n*i+j].im;
    }
  }

  // add spin-orbit contribution to the Hamiltonian
  // it has components in all quandrants of the martix
  for (i = 0; i < 4*n*n; i++) {
    hamiltonian[i].re += spinOrbitMatrix[i].re;
    hamiltonian[i].im += spinOrbitMatrix[i].im;
  }

  return;
}

/****************************************************************************/

void calcKineticEnergyMatrix(complexnumber *kineticEnergyMat, vector *basisStates, vector kVector, param params) {
  FILE *pf;
  int i, j, n = params.nBasisVectors; 
  double preFactor;
  vector kPlusG;

  preFactor = ( params.kineticEnergyScaling * sqr(HBAR) / (2.0 * MASS) );
 
  // The kinetic energy will only adds to the diagonal elements of the hamiltonian 
  // there is a delta function on both G (the basis vector) and the spin  
  for (i = 0; i < n; i++) {
     kPlusG = retAddedVectors(basisStates[i], kVector);
     kineticEnergyMat[2*n*i+i].re = sqr(kPlusG.mag) * preFactor;
     kineticEnergyMat[2*n*n+n+2*n*i+i].re = sqr(kPlusG.mag) * preFactor; 
  }    

  return;
}

/****************************************************************************/

void calcLocalPotentialEnergyMatrix(complexnumber *localPotential, vector *basisStates, atom *atoms, param params) { 
  FILE *pf;
  int i, j, k, s;
  double *gDiffDotTau;
  vector gDiff; 
  complexnumber *structFact; 

  int n = params.nBasisVectors;
  zeroComplexMatrix(localPotential, 2*n, 2*n);
  gDiffDotTau = (double *) calloc(params.nAtoms, sizeof(double));
  structFact = (complexnumber *) calloc(params.nAtoms, sizeof(complexnumber));

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      gDiff = retSubtractedVectors(basisStates[i], basisStates[j]);
      for (k = 0; k < params.nAtoms; k++) { 
        gDiffDotTau[k] = retDotProduct(gDiff, atoms[k].pos);
        structFact[k].re = (1.0 / (double) (params.nAtoms)) * cos(gDiffDotTau[k]);
        structFact[k].im = (-1.0 / (double) (params.nAtoms)) * sin(gDiffDotTau[k]);
        atoms[k].fF = calcPot(gDiff.mag, atoms[k].ppParams);
        // local potential energy has a delta function on spin -> block-diagonal     
        for (s = 0; s < 2; s++) {
           localPotential[2*n*n*s+s*n+2*n*i+j].re += atoms[k].fF * structFact[k].re;
           localPotential[2*n*n*s+s*n+2*n*i+j].im += atoms[k].fF * structFact[k].im; 
        }
      }  
    }
  }

  free(gDiffDotTau);
  free(structFact);
  
  return;
}

/****************************************************************************/

void calcSpinOrbitPotentialMatrix(complexnumber *spinOrbitMatrix, double *vSpinOrbit, vector *basisStates, vector kVector, int kindex, atom *atoms, param params) {
  int i, j, k, s;
  double *gDiffDotTau, preFactor, *lambda;
  vector ibsPlusK, jbsPlusK, gCrossProduct, gDiff;
  vector upUpRe, upUpIm, upDownRe, upDownIm, downUpRe, downUpIm, downDownRe, downDownIm;
  complexnumber gcpDotUpUp, gcpDotDownDown, gcpDotDownUp, gcpDotUpDown;
  complexnumber *structFact;

  gDiffDotTau = (double *) calloc(params.nAtoms, sizeof(double));
  structFact = (complexnumber *) calloc(params.nAtoms, sizeof(complexnumber));
  lambda = (double *) calloc(params.nAtoms, sizeof(double));

  // -i*<s|SpinMatrices|s'> vectors:
  upUpRe = retZeroVector();
  upUpIm.x = 0.0, upUpIm.y = 0.0, upUpIm.z = -0.5;
  downDownRe = retZeroVector();
  downDownIm.x = 0.0, downDownIm.y = 0.0, downDownIm.z = 0.5;
  upDownRe.x = 0.0, upDownRe.y = -0.5, upDownRe.z = 0.0;
  upDownIm.x = -0.5, upDownIm.y = 0.0, upDownIm.z = 0.0;
  downUpRe.x = 0.0, downUpRe.y = 0.5, downUpRe.z = 0.0;
  downUpIm.x = -0.5, downUpIm.y = 0.0, downUpIm.z = 0.0;

  int n = params.nBasisVectors;
  zeroComplexMatrix(spinOrbitMatrix, 2*n, 2*n); 
  
  for (i = 0; i < n; i++) { // i = k
    for (j = 0; j < n; j++) { // j = k'
      gDiff = retSubtractedVectors(basisStates[i], basisStates[j]);
      // TODO (Mar 22 2019): check if gDiff should be j-i instead of i-j
      ibsPlusK = retAddedVectors(basisStates[i], kVector);
      jbsPlusK = retAddedVectors(basisStates[j], kVector);    
      // avoid dividing by 0 in spinOrbitPrefactor
      if (ibsPlusK.mag < EPS || jbsPlusK.mag < EPS) {
        continue;
      }
      else {
        preFactor = 12.0 * PIE / (ibsPlusK.mag * jbsPlusK.mag);
        gCrossProduct = retCrossProduct(ibsPlusK, jbsPlusK);
        gcpDotUpUp.re = retDotProduct(upUpRe, gCrossProduct); // is always 0
        gcpDotUpUp.im = retDotProduct(upUpIm, gCrossProduct);
        gcpDotDownDown.re = retDotProduct(downDownRe, gCrossProduct); // is always 0
        gcpDotDownDown.im = retDotProduct(downDownIm, gCrossProduct);
        gcpDotUpDown.re = retDotProduct(upDownRe, gCrossProduct);
        gcpDotUpDown.im = retDotProduct(upDownIm, gCrossProduct);
        gcpDotDownUp.re = retDotProduct(downUpRe, gCrossProduct);
        gcpDotDownUp.im = retDotProduct(downUpIm, gCrossProduct); 
        for (k = 0; k < params.nAtoms; k++) { 
          gDiffDotTau[k] = retDotProduct(gDiff, atoms[k].pos);
          structFact[k].re = (1.0 / (double) (params.nAtoms)) * cos(gDiffDotTau[k]);
          structFact[k].im = (-1.0 / (double) (params.nAtoms)) * sin(gDiffDotTau[k]);
          lambda[k] = atoms[k].spinOrbit * vSpinOrbit[kindex*n*n+i*n+j];
          // up-up (top left) block of the matrix
          spinOrbitMatrix[2*n*i+j].re += preFactor*gcpDotUpUp.re*lambda[k]*structFact[k].re; // is always 0
          spinOrbitMatrix[2*n*i+j].im += preFactor*gcpDotUpUp.im*lambda[k]*structFact[k].im;
          // down-down (bottom right) block of the matrix
          spinOrbitMatrix[2*n*n+n+2*n*i+j].re += preFactor*gcpDotDownDown.re*lambda[k]*structFact[k].re; // is always 0
          spinOrbitMatrix[2*n*n+n+2*n*i+j].im += preFactor*gcpDotDownDown.im*lambda[k]*structFact[k].im;
          // up-down (top right) block of the matrix
          spinOrbitMatrix[n+2*n*i+j].re += preFactor*gcpDotUpDown.re*lambda[k]*structFact[k].re;
          spinOrbitMatrix[n+2*n*i+j].im += preFactor*gcpDotUpDown.im*lambda[k]*structFact[k].im;
          // down-up (bottom left) block of the matrix
          spinOrbitMatrix[2*n*n+2*n*i+j].re += preFactor*gcpDotDownUp.re*lambda[k]*structFact[k].re;
          spinOrbitMatrix[2*n*n+2*n*i+j].im += preFactor*gcpDotDownUp.im*lambda[k]*structFact[k].im;
        }
        // TODO: check if dot product in outside or inside of loop over the atoms matters
      }
    }
  }

  free(gDiffDotTau);
  free(structFact);
  free(lambda);

  return;
}

/****************************************************************************/
// Calculates Vso(K,K') = 1/cellVolume * integral from 0 t0 infinity of
// dr*r^2*j1(Kr)*exp^(-(r/0.7)^2)*j1(K'r) where j1 is the 1st bessel function
// K = kpoint + basisVector and exp^(-(r/0.7)^2) is the  spin-orbit potential 
// (Hybertsen & Louie, PRB, 34, 2920 (1986))
// The nKPoint*nBasisVectors^2 of these are stored in vSpinOrbit

void calcVSpinOrbit(double *vSpinOrbit, vector *basisStates, vector *kPoint, param params) {
  int i, j, k, gp;
  double r, eps, sum;
  int n = params.nBasisVectors;
  double dr = TWOPI/(100.0*basisStates[n-1].mag);
  double width = 0.7;
  double rCut = sqrt(width*width*16.0*log(10.0));
  int nCut = (int)(rCut/dr);
  vector gikp, gjkp;

  eps = dr / 10.0;

  for (k = 0; k < params.nKPoints; k++) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        gikp = retAddedVectors(basisStates[i], kPoint[k]);
        gjkp = retAddedVectors(basisStates[j], kPoint[k]);
        sum = 0.0;
        // perform numerical integration 
        for (gp = 1; gp < nCut; gp++) {
          r = ((double)(gp) * dr);
          sum += (sqr(r) * calcBessel(gikp.mag*r, eps) * exp(-(sqr(r/width))) * calcBessel(gjkp.mag*r, eps) * dr);                       
        }
        vSpinOrbit[k*n*n+i*n+j] = sum / (params.cellVolume);
      }
    }
  }

  return;
}

/****************************************************************************/

double calcPot(double q, double *param) {
  double pot = (param[0]*(q*q - param[1]) / (param[2] * exp(param[3]*q*q) + 1.0));
  return pot;
}

/****************************************************************************/

void zeroComplexMatrix(complexnumber *mat, int n1, int n2) {
  int i, n = n1*n2;

  for (i = 0; i < n; i++) {
    mat[i].re = 0.0;
    mat[i].im = 0.0;
  }

  return;
}

/****************************************************************************/
