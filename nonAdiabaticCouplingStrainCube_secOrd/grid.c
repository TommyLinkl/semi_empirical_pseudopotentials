/*****************************************************************************/

#include "qp.h"

/****************************************************************************/
/* initialized the real-space grid */
/****************************************************************************/
long initRSpaceGrid(grid3d *rSpaceGrid, gridPoint3d *rSpaceGP, lParams lPar) {

    FILE *pf;
    long iGrid, iX, iY, iZ, nAtoms;
    double x, y, z;

    // r-space grid point will be indexed by 0
    rSpaceGrid->index = 0;
    rSpaceGrid->gP = rSpaceGP;

    // determine the dimensions of the r-space grid required 
    pf = fopen("conf.par" , "r");
    fscanf(pf, "%ld", &nAtoms);
    determineRSpaceGridSize(rSpaceGrid, nAtoms, pf);
    fclose(pf);

    // fill in the grid point structures for all grid points
    // the grid points go in the order of z->y->x (i.e. x changes the quickest)
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
/* determines size of grid and spacing between grid points from conf file */
/****************************************************************************/
long determineRSpaceGridSize(grid3d *rSpaceGrid, long nAtoms, FILE *confFilePointer) {
    long iAtom;
    double rx[nAtoms], ry[nAtoms], rz[nAtoms];
    atm_st atom[nAtoms];
    
    // Read in the atom positions to be used in determining the size of the r-space grid
    readConfFile(rx, ry, rz, atom, nAtoms, confFilePointer);
    
    // box dimensions
    rSpaceGrid->maxPos.x = rint(0.5 * get_dot_ligand_size_z(rx, nAtoms) + 5.0);
    rSpaceGrid->maxPos.y = rint(0.5 * get_dot_ligand_size_z(ry, nAtoms) + 5.0);
    rSpaceGrid->maxPos.z = rint(0.5 * get_dot_ligand_size_z(rz, nAtoms) + 5.0);
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
