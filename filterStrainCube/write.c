/****************************************************************************/
//
// This file does the printing for the program 
//
/****************************************************************************/

#include "fd.h"

/****************************************************************************/
// prints out current time to stdout 
 
void writeCurrentTime(FILE *pf) {
  time_t startTime;

  startTime = time(NULL);
  fprintf(pf, "%s\n", ctime(&startTime));

  return;
}

/****************************************************************************/

void writeSeparation(FILE *pf) {
  fprintf(pf, "\n******************************************************************************\n\n");
  
  return;
}

/****************************************************************************/

void writeCubeFile(double *rho, double xmin, double ymin, double zmin, 
        double dx, double dy, double dz, long nx, long ny, long nz) {

    FILE *pf, *pConfFile;
    long iGrid, iX, iY, iZ, iYZ, iXY, nAtoms, atomType;
    double x, y, z;
    char line[80], atomSymbol[10];

    pConfFile = fopen("conf.dat", "r"); //conf.par
    fscanf(pConfFile, "%ld", &nAtoms);
    pf = fopen("localPot.cube", "w");
    fprintf(pf, "CUBE FILE\n");
    fprintf(pf, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");

    fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nAtoms, xmin, ymin, zmin);
    fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nx, dx, 0.0, 0.0);
    fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ny, 0.0, dy, 0.0);
    fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nz, 0.0, 0.0, dz);
    
    fgets(line, 80, pConfFile); 
    while(fgets(line, 80, pConfFile) != NULL) {
        sscanf(line, "%2s %lf %lf %lf %ld %*lf", atomSymbol, &x, &y, &z, &atomType);
        if (! strcmp(atomSymbol, "Cd")) { 
            atomType = 48;
        }
        else if (! strcmp(atomSymbol, "S")) {  
            atomType = 16;
        }
        else if (! strcmp(atomSymbol, "Se")) { 
            atomType = 34;
        }
        else if (! strcmp(atomSymbol, "Zn")) {
            atomType = 30;
        }
        else if (! strcmp(atomSymbol, "Te")) {
            atomType = 52;
        }
        else if (! strcmp(atomSymbol, "C")) {
            atomType = 6;
        }
        else if (! strcmp(atomSymbol, "Si")) {
            atomType = 14;
        }
        else if (! strcmp(atomSymbol, "Cs")) {
            atomType = 55;
        }
        else if (! strcmp(atomSymbol, "Pb")) {
            atomType = 82;
        }
        else if (! strcmp(atomSymbol, "I")) {
            atomType = 53;
        }
        else if (! strcmp(atomSymbol, "In")) {
            atomType = 49;
        }
        else if (! strcmp(atomSymbol, "As")) {
            atomType = 33;
        }
        else if (! strcmp(atomSymbol, "Ga")) {
            atomType = 31;
        }
        else if (! strcmp(atomSymbol, "P")) {
            atomType = 15;
        }
        else { 
            atomType = 1; 
        }
        fprintf(pf, "%5li%12.6f%12.6f%12.6f%12.6f\n", atomType, 0.0, x, y, z);
    }
    
    for (iX = 0; iX < nx; iX++) {
        for (iY = 0; iY < ny; iY++) {
            for (iZ = 0; iZ < nz; iZ++) {
                iGrid = iZ*ny*nx + iY*nx + iX;
                fprintf(pf, "%g ", rho[iGrid]);
                if (iZ % 6 == 5) {fprintf(pf, "\n"); }
            }
            fprintf(pf, "\n");
        }
    }

    fclose(pConfFile);
    fclose(pf);

    return;
}
