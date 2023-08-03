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

void writeCubeFile(double *rho, par_st par, long_st ist, char *fileName) {
  FILE *pf, *pConfFile,*pConfFile2;
  long iGrid, iX, iY, iZ, iYZ, nAtoms, atomType;
  double x, y, z;
  char line[80], atomSymbol[10];


  pConfFile = fopen("conf.dat", "r");
  pConfFile2 = fopen("conf.dat", "r"); // has to be .par like
  fscanf(pConfFile2, "%ld", &nAtoms);
  pf = fopen(fileName, "w");
  fprintf(pf, "CUBE FILE\n");
  fprintf(pf, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Y\n");
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nAtoms, par.xmin, par.ymin, par.zmin);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ist.nx, par.dx, 0.0, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ist.ny, 0.0, par.dy, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ist.nz, 0.0, 0.0, par.dz);
  fgets(line, 80, pConfFile); 
  while(fgets(line, 80, pConfFile) != NULL) {
    sscanf(line, "%2s %lf %lf %lf %ld %*lf", atomSymbol, &x, &y, &z, &atomType);
    
    //TODO: make sane atom handling.... 
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
  for (iX = 0; iX < ist.nx; iX++) {
    for (iY = 0; iY < ist.ny; iY++) {
        for (iZ = 0; iZ < ist.nz; iZ++) {
        iYZ = ist.nx * (ist.ny * iZ + iY);
        iGrid = iYZ + iX;
        fprintf(pf, "%g ", rho[iGrid]);
        if (iZ % 6 == 5) {
          fprintf(pf, "\n");
        }
      }
      fprintf(pf, "\n");
    }
  }
  fclose(pConfFile);
  fclose(pf);

  return;
}


void writeCubeFile_cubicUnitCell(double *rho, par_st par, long_st ist, char *fileName, double box_xmin, double box_ymin, double box_zmin, double box_length) {
  FILE *pf, *pConfFile,*pConfFile2;
  long iGrid, iX, iY, iZ, iYZ, nAtoms, atomType;
  double grid_x, grid_y, grid_z, dist_r; 
  double x, y, z;
  char line[80], atomSymbol[10];
  double box_xmax, box_ymax, box_zmax; 

  box_xmax = box_xmin + box_length; 
  box_ymax = box_ymin + box_length; 
  box_zmax = box_zmin + box_length;

  pConfFile = fopen("conf.dat", "r");
  pConfFile2 = fopen("conf.dat", "r"); // has to be .par like
  fscanf(pConfFile2, "%ld", &nAtoms);
  pf = fopen(fileName, "w");
  fprintf(pf, "CUBE FILE\n");
  fprintf(pf, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nAtoms, par.xmin, par.ymin, par.zmin);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ist.nx, par.dx, 0.0, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ist.ny, 0.0, par.dy, 0.0);
  fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ist.nz, 0.0, 0.0, par.dz);
  fgets(line, 80, pConfFile); 
  while(fgets(line, 80, pConfFile) != NULL) {
    sscanf(line, "%2s %lf %lf %lf %ld %*lf", atomSymbol, &x, &y, &z, &atomType);
    
    //TODO: make sane atom handling.... 
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
    if ((x>=box_xmin) && (x<=box_xmax) && (y>=box_ymin) && (y<=box_ymax) && (z>=box_zmin) && (z<=box_zmax)) {
      fprintf(pf, "%5li%12.6f%12.6f%12.6f%12.6f\n", atomType, 0.0, x, y, z);
    }
  }
  for (iX = 0; iX < ist.nx; iX++) {
    for (iY = 0; iY < ist.ny; iY++) {
      for (iZ = 0; iZ < ist.nz; iZ++) {
        iYZ = ist.nx * (ist.ny * iZ + iY);
        iGrid = iYZ + iX;
        if ((iX<=ist.nx/2-12) || (iX>=ist.nx/2+12) || (iY<=ist.ny/2-12) || (iY>=ist.ny/2+12) || (iZ<=ist.nz/2-12) || (iZ>=ist.nz/2+12)) {fprintf(pf, "%g ", 0.0);} //8
        else {fprintf(pf, "%g ", rho[iGrid]);}
        if (iZ % 6 == 5) {fprintf(pf, "\n");}
      }
      fprintf(pf, "\n");
    }
  }
  fclose(pConfFile);
  fclose(pf);

  return;
}



/****************************************************************************/
// Dipti version
// void writeCubeFile(double *rho, double xmin, double ymin, double zmin, 
//         double dx, double dy, double dz, long nx, long ny, long nz) {

//     FILE *pf, *pConfFile;
//     long iGrid, iX, iY, iZ, iYZ, nAtoms, atomType;
//     double x, y, z;
//     char line[80], atomSymbol[10];

//     pConfFile = fopen("conf.par", "r");
//     fscanf(pConfFile, "%ld", &nAtoms);
//     pf = fopen("localPot.cube", "w");
//     fprintf(pf, "CUBE FILE\n");
//     fprintf(pf, "OUTER LOOP: Z, MIDDLE LOOP: Y, INNER LOOP: X\n");

//     fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nAtoms, xmin, ymin, zmin);
//     fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nz, 0.0, 0.0, dz);
//     fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ny, 0.0, dy, 0.0);
//     fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nx, dx, 0.0, 0.0);
//     fgets(line, 80, pConfFile); 
//     while(fgets(line, 80, pConfFile) != NULL) {
//         sscanf(line, "%2s %lf %lf %lf %d %*lf", &atomSymbol, &x, &y, &z, &atomType);
//         if (! strcmp(atomSymbol, "Cd")) { 
//             atomType = 48;
//         }
//         else if (! strcmp(atomSymbol, "S")) {  
//             atomType = 16;
//         }
//         else if (! strcmp(atomSymbol, "Se")) { 
//             atomType = 34;
//         }
//         else if (! strcmp(atomSymbol, "Zn")) {
//             atomType = 30;
//         }
//         else if (! strcmp(atomSymbol, "Te")) {
//             atomType = 52;
//         }
//         else if (! strcmp(atomSymbol, "C")) {
//             atomType = 6;
//         }
//         else if (! strcmp(atomSymbol, "Si")) {
//             atomType = 14;
//         }
//         else if (! strcmp(atomSymbol, "Cs")) {
//             atomType = 55;
//         }
//         else if (! strcmp(atomSymbol, "Pb")) {
//             atomType = 82;
//         }
//         else if (! strcmp(atomSymbol, "I")) {
//             atomType = 53;
//         }
//         else { 
//             atomType = 1; 
//         }
//         fprintf(pf, "%5i%12.6f%12.6f%12.6f%12.6f\n", atomType, 0.0, x, y, z);
//     }
    
//     for (iZ = 0; iZ < nz; iZ++) {
//         for (iY = 0; iY < ny; iY++) {
//             iYZ = nx * (ny * iZ + iY);
//             for (iX = 0; iX < nx; iX++) {
//                 iGrid = iYZ + iX;
//                 fprintf(pf, "%g ", rho[iGrid]);
//                 if (iX % 6 == 5) {
//                     fprintf(pf, "\n");
//                 }
//             }
//             fprintf(pf, "\n");
//         }
//     }
//     fclose(pConfFile);
//     fclose(pf);

//     return;
// }
