#include "qp.h"
#include <omp.h>

/******************************************************************************/
/* calculate derivative couplings (i.e. el-ph coupling for adiabatic states)
 * between quasiparticle states
 * int(psi_r * dPot/dR * psi_s) / (E_s - E_r) for all atoms */
/******************************************************************************/
/*
void calcDerivativeCouplings(vector *Vab, vector *Vij, vector *dPotdR, atom *atoms, 
        grid1d dAtomicPPGrid, nonintQP *holeQP, nonintQP *elecQP, 
        grid3d rSpaceGrid, lParams lPar) {

    FILE *pf;
    long i, iR, iS, iGrid, iAtom, atomType, nTotalPPValues;
    double *dPotValues, *dPotGridPointValues;
    double energyR, energyS, rho, dVOverE;
    double *psiR, *psiS;
    vector sum, tmpVector;

    // dynamically allocate memory
    nTotalPPValues = dAtomicPPGrid.nGridPoints*lPar.nSCAtomTypes;
    if ((dPotValues = (double *) calloc(nTotalPPValues, sizeof(double))) == NULL) memoryError("dPotValues");
    if ((dPotGridPointValues = (double *) calloc(nTotalPPValues, sizeof(double))) == NULL) memoryError("dPotGridPointValues");

    // get potential grid points and values
    for (iGrid = 0; iGrid < nTotalPPValues; iGrid++) {
        dPotGridPointValues[iGrid] = dAtomicPPGrid.gP[iGrid].pos;
        dPotValues[iGrid] = dAtomicPPGrid.gP[iGrid].localVal;
    }

    // compute nonadiabatic matrix elements for each atom
    for (iAtom = 0; iAtom < lPar.nAtoms; iAtom++) {
        if (isAPassivationSymbol(atoms[iAtom].symbol)) {
            continue;
        }
        else {
            atomType = atoms[iAtom].type;
            i = atomType*dAtomicPPGrid.nGridPoints;

            calcdPotentialdR(dPotdR, rSpaceGrid, atoms[iAtom].pos, &(dPotValues[i]), &(dPotGridPointValues[i]),
                             dAtomicPPGrid.stepSize, dAtomicPPGrid.nGridPoints);

            // calculate electron matrix elements
            for (iR = 0; iR < lPar.nElecs; iR++) {
                // diagonal elements are 0
                Vab[iAtom*lPar.nElecs*lPar.nElecs + iR*lPar.nElecs + iR] = retZeroVector();

                energyR = elecQP[iR].energy;
                psiR = elecQP[iR].psi;
                for (iS = (iR+1); iS < lPar.nElecs; iS++) {
                    energyS = elecQP[iS].energy;
                    psiS = elecQP[iS].psi;
                    dVOverE = rSpaceGrid.dV / (energyS - energyR);

                    // perform integral
                    sum = retZeroVector();
                    for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
                        rho = psiR[iGrid]*psiS[iGrid];
                        tmpVector = retScaledVector(dPotdR[iGrid], rho);
                        sum = retAddedVectors(sum, tmpVector);
                    }
                    Vab[iAtom*lPar.nElecs*lPar.nElecs + iR*lPar.nElecs + iS] = retScaledVector(sum, dVOverE);
                    // Vba = -Vab
                    Vab[iAtom*lPar.nElecs*lPar.nElecs + iS*lPar.nElecs + iR] = retScaledVector(sum, -dVOverE);
                }
            }

            // calculate hole matrix elements
            for (iR = 0; iR < lPar.nHoles; iR++) {
                // diagonal elements are 0
                Vij[iAtom*lPar.nHoles*lPar.nHoles + iR*lPar.nHoles + iR] = retZeroVector();

                energyR = holeQP[iR].energy;
                psiR = holeQP[iR].psi;
                for (iS = (iR+1); iS < lPar.nHoles; iS++) {
                    energyS = holeQP[iS].energy;
                    psiS = holeQP[iS].psi;
                    dVOverE = rSpaceGrid.dV / (energyS - energyR);

                    // perform integral
                    sum = retZeroVector();
                    for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
                        rho = psiR[iGrid]*psiS[iGrid];
                        tmpVector = retScaledVector(dPotdR[iGrid], rho);
                        sum = retAddedVectors(sum, tmpVector);
                    }
                    Vij[iAtom*lPar.nHoles*lPar.nHoles + iR*lPar.nHoles + iS] = retScaledVector(sum, dVOverE);
                    // Vji = -Vij
                    Vij[iAtom*lPar.nHoles*lPar.nHoles + iS*lPar.nHoles + iR] = retScaledVector(sum, -dVOverE);
                }
            }

        }
    }

    if (lPar.nElecs > 0) {
        // write matrix elements to file
        pf = fopen("Vab.dat", "w");
        for (iAtom = 0; iAtom < lPar.nSCAtoms; iAtom++) {
            for (iR = 0; iR < lPar.nElecs; iR++) {
                for (iS = 0; iS < lPar.nElecs; iS++) {
                    i = iAtom*lPar.nElecs*lPar.nElecs + iR*lPar.nElecs + iS;
                    fprintf(pf, "%ld %s %lg %lg %lg %ld % .8f %ld % .8f %lg %lg %lg %lg\n", 
                            iAtom, atoms[iAtom].symbol, atoms[iAtom].pos.x, 
                            atoms[iAtom].pos.y, atoms[iAtom].pos.z, 
                            iR, elecQP[iR].energy, iS, elecQP[iS].energy, 
                            Vab[i].x, Vab[i].y, Vab[i].z, Vab[i].mag);
                }
            }
        }
        fclose(pf);
    }

    if (lPar.nHoles > 0) {
        pf = fopen("Vij.dat", "w");
        for (iAtom = 0; iAtom < lPar.nSCAtoms; iAtom++) {
            for (iR = 0; iR < lPar.nHoles; iR++) {
                for (iS = 0; iS < lPar.nHoles; iS++) {
                    i = iAtom*lPar.nHoles*lPar.nHoles + iR*lPar.nHoles + iS;
                    fprintf(pf, "%ld %s %lg %lg %lg %ld % .8f %ld % .8f %lg %lg %lg %lg\n", 
                            iAtom, atoms[iAtom].symbol, atoms[iAtom].pos.x, 
                            atoms[iAtom].pos.y, atoms[iAtom].pos.z, 
                            iR, holeQP[iR].energy, iS, holeQP[iS].energy, 
                            Vij[i].x, Vij[i].y, Vij[i].z, Vij[i].mag);
                }
            }
        }
        fclose(pf);
    }

    // free dynamically allocated memory
    free(dPotValues); free(dPotGridPointValues);

    return;
}
*/

/******************************************************************************/
/* calculate electron-phonon matrix elements (i.e. el-ph coupling for diabatic 
 * states) between quasiparticle states
 * int(psi_r * dPot/dR * psi_s) for all atoms 
 * computes diagonal matrix elements (for reorganization energy / HR parameter) */
/******************************************************************************/
void calcElPh(vector *Vab, vector *Vij, vector *dPotdR, atom *atoms, atom *atomNeighbors,
        grid1d atomicPPGrid, grid1d dAtomicPPGrid, nonintQP *holeQP, nonintQP *elecQP, 
        grid3d rSpaceGrid, double *strainScale, vector *strainScaleDeriv, lParams lPar) {

    FILE *pf;
    long i, j, iR, iS, iGrid, iAtom, iNeighbor, neighborIndex, atomType, nTotalPPValues;
    double *potValues, *potGridPointValues;
    double *dPotValues, *dPotGridPointValues;
    double rho;
    double *psiR, *psiS;
    vector sum, tmpVector;
    
    printf("In calcElPh, starting allocating memory.\n"); 
    // dynamically allocate memory
    nTotalPPValues = dAtomicPPGrid.nGridPoints*lPar.nSCAtomTypes;
    if ((potValues = (double *) calloc(nTotalPPValues, sizeof(double))) == NULL) memoryError("potValues");
    if ((potGridPointValues = (double *) calloc(nTotalPPValues, sizeof(double))) == NULL) memoryError("potGridPointValues");
    if ((dPotValues = (double *) calloc(nTotalPPValues, sizeof(double))) == NULL) memoryError("dPotValues");
    if ((dPotGridPointValues = (double *) calloc(nTotalPPValues, sizeof(double))) == NULL) memoryError("dPotGridPointValues");
    printf("In calcElPh, finished with allocating memory. \n"); 

    // get potential grid points and values
    for (iGrid = 0; iGrid < nTotalPPValues; iGrid++) {
        potGridPointValues[iGrid] = atomicPPGrid.gP[iGrid].pos;
        potValues[iGrid] = atomicPPGrid.gP[iGrid].localVal;
    }
    // get potential derivative grid points and values
    for (iGrid = 0; iGrid < nTotalPPValues; iGrid++) {
        dPotGridPointValues[iGrid] = dAtomicPPGrid.gP[iGrid].pos;
        dPotValues[iGrid] = dAtomicPPGrid.gP[iGrid].localVal;
    }

    // compute nonadiabatic matrix elements for each atom
    for (iAtom = 0; iAtom < lPar.nAtoms; iAtom++) {
        printf("Atom number: %d \n", iAtom);
        if (isAPassivationSymbol(atoms[iAtom].symbol)) {
            continue;
        }
        else {

            atomType = atoms[iAtom].type;
            i = atomType*dAtomicPPGrid.nGridPoints;

            // zero dPotdR
            for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
                dPotdR[iGrid] = retZeroVector();
            }
            // derivative of iAtom's pseudopotential wrt iAtom's position (on grid)
            calcdPotentialdR(dPotdR, rSpaceGrid, atoms[iAtom].pos, &(dPotValues[i]), &(dPotGridPointValues[i]),
                             dAtomicPPGrid.stepSize, dAtomicPPGrid.nGridPoints, strainScale[iAtom]);
            // // iAtom's neighbor's pseudopotentials (on grid)
            for (iNeighbor = 0; iNeighbor < 4; iNeighbor++) {
                neighborIndex = atomNeighbors[4*iAtom+iNeighbor].index;
                if (! isAPassivationSymbol(atoms[neighborIndex].symbol)) {
                    j = atoms[neighborIndex].type*atomicPPGrid.nGridPoints;
                    calcPotential(dPotdR, rSpaceGrid, atoms[neighborIndex].pos, &(potValues[j]), &(potGridPointValues[j]),
                            atomicPPGrid.stepSize, atomicPPGrid.nGridPoints, strainScaleDeriv[4*iAtom+iNeighbor]);
                }
            }

            printf("Starting to compute electron matrix elements\n");
            // calculate electron matrix elements
            for (iR = 0; iR < lPar.nElecs; iR++) {
                psiR = elecQP[iR].psi;

                // diagaonal element
                sum = retZeroVector();
                for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
                    rho = psiR[iGrid]*psiR[iGrid];
                    tmpVector = retScaledVector(dPotdR[iGrid], rho);
                    sum = retAddedVectors(sum, tmpVector);
                }
                Vab[iAtom*lPar.nElecs*lPar.nElecs + iR*lPar.nElecs + iR] = retScaledVector(sum, rSpaceGrid.dV);

                #pragma omp parallel for private(psiS, sum, iGrid, rho, tmpVector)
                // non-diagonal elements
                for (iS = (iR+1); iS < lPar.nElecs; iS++) {
                    int tid, nthreads;
                    nthreads = omp_get_num_threads();
                    // printf("Total threads: %d\n", nthreads);
                    tid = omp_get_thread_num();
                    // printf("Matrix element %ld %ld for atom %ld computed on thread %d\n", iR, iS, iAtom, tid);
                   
                    psiS = elecQP[iS].psi;

                    sum = retZeroVector();
                    for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
                        rho = psiR[iGrid]*psiS[iGrid];
                        tmpVector = retScaledVector(dPotdR[iGrid], rho);
                        sum = retAddedVectors(sum, tmpVector);
                    }
                    Vab[iAtom*lPar.nElecs*lPar.nElecs + iR*lPar.nElecs + iS] = retScaledVector(sum, rSpaceGrid.dV);
                    // Vba = Vab
                    Vab[iAtom*lPar.nElecs*lPar.nElecs + iS*lPar.nElecs + iR] = retScaledVector(sum, rSpaceGrid.dV);
                }
            }

            printf("Starting to compute hole matrix elements\n");
            // calculate hole matrix elements
            for (iR = 0; iR < lPar.nHoles; iR++) {
                psiR = holeQP[iR].psi;

                // diagonal element
                sum = retZeroVector();
                for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
                    rho = psiR[iGrid]*psiR[iGrid];
                    tmpVector = retScaledVector(dPotdR[iGrid], rho);
                    sum = retAddedVectors(sum, tmpVector);
                }
                Vij[iAtom*lPar.nHoles*lPar.nHoles + iR*lPar.nHoles + iR] = retScaledVector(sum, rSpaceGrid.dV);

                #pragma omp parallel for private(psiS, sum, iGrid, rho, tmpVector)
                // non-diagonal elements
                for (iS = (iR+1); iS < lPar.nHoles; iS++) {
                    int tid, nthreads;
                    nthreads = omp_get_num_threads();
                    // printf("Total threads: %d\n", nthreads);
                    tid = omp_get_thread_num();
                    // printf("Matrix element %ld %ld for atom %ld computed on thread %d\n", iR, iS, iAtom, tid);

                    psiS = holeQP[iS].psi;

                    sum = retZeroVector();
                    for (iGrid = 0; iGrid < rSpaceGrid.nGridPoints; iGrid++) {
                        rho = psiR[iGrid]*psiS[iGrid];
                        tmpVector = retScaledVector(dPotdR[iGrid], rho);
                        sum = retAddedVectors(sum, tmpVector);
                    }
                    Vij[iAtom*lPar.nHoles*lPar.nHoles + iR*lPar.nHoles + iS] = retScaledVector(sum, rSpaceGrid.dV);
                    // Vji = Vij
                    Vij[iAtom*lPar.nHoles*lPar.nHoles + iS*lPar.nHoles + iR] = retScaledVector(sum, rSpaceGrid.dV);
                }
            }

        }
    }
    printf("Done compute nonadiabatic matrix elements for all atoms \n");

    printf("Starting to write nonadiabatic matrix elements to files \n");
    if (lPar.nElecs > 0) {
        // write matrix elements to file
        pf = fopen("Vab-diabatic.dat", "w");
        for (iAtom = 0; iAtom < lPar.nSCAtoms; iAtom++) {
            for (iR = 0; iR < lPar.nElecs; iR++) {
                for (iS = 0; iS < lPar.nElecs; iS++) {
                    i = iAtom*lPar.nElecs*lPar.nElecs + iR*lPar.nElecs + iS;
                    fprintf(pf, "%ld %s %lg %lg %lg %lg ",
                        iAtom, atoms[iAtom].symbol, atoms[iAtom].pos.x, 
                        atoms[iAtom].pos.y, atoms[iAtom].pos.z, atoms[iAtom].pos.mag);
                    fprintf(pf, "%ld % .8f %ld % .8f ",
                        iR, elecQP[iR].energy, iS, elecQP[iS].energy);
                    fprintf(pf, "%lg %lg %lg %lg\n",
                        Vab[i].x, Vab[i].y, Vab[i].z, Vab[i].mag);
                }
            }
        }
        fclose(pf);
    }

    if (lPar.nHoles > 0) {
        pf = fopen("Vij-diabatic.dat", "w");
        for (iAtom = 0; iAtom < lPar.nSCAtoms; iAtom++) {
            for (iR = 0; iR < lPar.nHoles; iR++) {
                for (iS = 0; iS < lPar.nHoles; iS++) {
                    i = iAtom*lPar.nHoles*lPar.nHoles + iR*lPar.nHoles + iS;
                    fprintf(pf, "%ld %s %lg %lg %lg %lg ", 
                        iAtom, atoms[iAtom].symbol, atoms[iAtom].pos.x, 
                        atoms[iAtom].pos.y, atoms[iAtom].pos.z, atoms[iAtom].pos.mag);
                    fprintf(pf, "%ld % .8f %ld % .8f ", 
                        iR, holeQP[iR].energy, iS, holeQP[iS].energy);
                    fprintf(pf, "%lg %lg %lg %lg\n", 
                        Vij[i].x, Vij[i].y, Vij[i].z, Vij[i].mag);
                }
            }
        }
        fclose(pf);
    }

    // free dynamically allocated memory
    free(potValues); free(potGridPointValues);
    free(dPotValues); free(dPotGridPointValues);

    return;
}
