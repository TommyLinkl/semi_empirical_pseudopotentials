#include "qp.h"

/******************************************************************************/
/* calculate derivative of pseudopotential wrt |r-R| */
/******************************************************************************/
void calcPseudopotentialDerivative(grid1d *dAtomicPPGrid, gridPoint1d *dAtomicPPGP,
                                   grid1d atomicPPGrid, long nAtomTypes) {

    FILE *pf;
    long i, iAtomType, iGrid;
    long nGridPoints;
    double dr;
    char fileName[100], atype[3];

    // real-space grid for derivative of pseudopotential indexed by 2
    dAtomicPPGrid->index = 2;
    dAtomicPPGrid->gP = dAtomicPPGP;
    dAtomicPPGrid->nGridPoints = atomicPPGrid.nGridPoints;
    dAtomicPPGrid->stepSize = atomicPPGrid.stepSize;
    nGridPoints = atomicPPGrid.nGridPoints;
    dr = atomicPPGrid.stepSize;

    // compute derivative of each pseudopotential
    for (iAtomType = 0; iAtomType < nAtomTypes; iAtomType++) {
        for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
            i = iAtomType*nGridPoints + iGrid;
            dAtomicPPGP[i].index = iGrid;
            dAtomicPPGP[i].pos = atomicPPGrid.gP[i].pos;

            // endpoints
            if (iGrid == 0 || iGrid == nGridPoints-1) {
                dAtomicPPGP[i].localVal = 0.;
            }
            // compute derivative via finite difference
            else {
                dAtomicPPGP[i].localVal = (atomicPPGrid.gP[i+1].localVal - atomicPPGrid.gP[i-1].localVal) / (2.*dr);
            }

        }
    }

    for (iAtomType = 0; iAtomType < nAtomTypes; iAtomType++) {

        strcpy(fileName, "dpot");
        assign_atom_type(atype, iAtomType);
        strcat(fileName, atype);
        strcat(fileName, ".dat");
        
        pf = fopen(fileName, "w");
        for (iGrid = 0; iGrid < nGridPoints; iGrid++) {
            i = iAtomType*nGridPoints + iGrid;
            fprintf(pf, "%.8f %.8f\n", dAtomicPPGP[i].pos, dAtomicPPGP[i].localVal);
        }
        fclose(pf);
    }

    return;
}

/******************************************************************************/
/* calculate atomic potential for all real-space grid points */
/******************************************************************************/
void calcPotential(vector *pot, grid3d rSpaceGrid, vector atomPosition, 
        double *potentialValues, double *potentialGridPointValues,
        double potGridStepSize, long nPotValues, vector strainScaleDeriv) {

    long iGrid;
    double potVal;
    double distR;
    vector delR;

    for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
        // grid point mimus atom (r-R)
        delR = retSubtractedVectors(rSpaceGrid.gP[iGrid].pos, atomPosition);
        // distance between grid point and atom (|r-R|)
        distR = delR.mag;

        // real-space cut-off -- if grid space is far from atom
        // potential is 0 at that grid point
        if (distR > 26.) {
            pot[iGrid] = retAddedVectors(pot[iGrid], retZeroVector());
        }
        else {
            // interpolate to find pseudopotential at distR (PP(|r-R|))
            potVal = interpolate(distR, potGridStepSize, potentialValues,
                               potentialGridPointValues, nPotValues);
            // multiply by derivative of strain scale
            pot[iGrid] = retAddedVectors(pot[iGrid], retScaledVector(strainScaleDeriv, potVal));
        }
    }

    return;
}

/******************************************************************************/
/* calculate derivative of potential wrt atomic coordinate (dPot/dR) 
 * calculates for all real-space grid points for a given atom */
/******************************************************************************/
void calcdPotentialdR(vector *dPotdR, grid3d rSpaceGrid, vector atomPosition,
                      double *potentialValues, double *potentialGridPointValues,
                      double potGridStepSize, long nPotValues, double strainScale) {

    long iGrid;
    double dPot;
    double distR;
    vector delR;

    for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
        // grid point mimus atom (r-R)
        delR = retSubtractedVectors(rSpaceGrid.gP[iGrid].pos, atomPosition);
        // distance between grid point and atom (|r-R|)
        distR = delR.mag;

        // real-space cut-off -- if grid space is far from atom
        // potential is 0 at that grid point
        if (distR > 26.) {
            dPotdR[iGrid] = retAddedVectors(dPotdR[iGrid], retZeroVector());
        }
        else {
            // interpolate to find derivative of pseudopotential at distR (dPP/d(|r-R|))
            dPot = interpolate(distR, potGridStepSize, potentialValues,
                               potentialGridPointValues, nPotValues);
            // multiply by d(|r-R|)/dR to get (dPot/dR) 
            // also multiply by strain scale 
            dPotdR[iGrid] = retAddedVectors(dPotdR[iGrid], retScaledVector(delR, -strainScale * dPot / (distR + EPS)));
        }
    }

    return;
}

/******************************************************************************/
/* linear interpolation for (derivative of) pseudopotential */
/******************************************************************************/
double interpolate(double evaluationPoint, double stepSize, double *potentialValues,
                   double *potentialGrid, long nPotValues) {

    long leftGridPointIdx;
    double slope, potAtLeftGridPoint, leftGridPoint, potAtEvalPoint;

    leftGridPointIdx = (long)(evaluationPoint/stepSize);
    leftGridPoint = potentialGrid[leftGridPointIdx];

    if (leftGridPointIdx > nPotValues-2) {
        potAtEvalPoint = 0.0;
    }
    else {
        potAtLeftGridPoint = potentialValues[leftGridPointIdx];
        slope = (potentialValues[leftGridPointIdx+1] - potAtLeftGridPoint) / stepSize;
        potAtEvalPoint = slope*(evaluationPoint - leftGridPoint) + potAtLeftGridPoint;
    }

    return potAtEvalPoint;
}
/******************************************************************************/
