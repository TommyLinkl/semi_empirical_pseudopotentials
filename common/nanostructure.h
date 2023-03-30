/****************************************************************************/
//
//
//
/****************************************************************************/

#ifndef NC_H
#define NC_H

/****************************************************************************/
/* These are the library functions that are used */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include "vector.h"
#include "atom.h"

/****************************************************************************/
/* Macro definitions: common unit conversion and multiplication schemes
  that help make the program more readable */

#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))
#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define AUTONM    0.05291772108
#define AUTOANG   0.5291772108
#define ANGTOAU   1.889725988579
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define HBAR      1.0
#define NDIM      3
#define MASS      1.0
#define EPS       1e-10

// All lattice constants are in Angstroms
#define ACdSeWZ 4.2989773 // Cd-Se bond length = 2.81 A, Eg = 1.714 eV
#define ACdSeZB 6.0583000 // Cd-Se bond length = 2.81 A, Eg =  
#define ACdSWZ  4.1302281 // Cd-S  bond length = 2.81 A, Eg = 2.50  eV
#define ACdSZB  5.8179695 // Cd-S  bond length = 2.81 A, Eg = 2.50  eV
#define ACdTeWZ 4.60      // TODO: find better value! Cd-Te bond length =  A, Eg =  eV 
#define ACdTeZB 6.48      // Cd-Te bond length = 2.81 A, Eg = 1.474 eV 
// www.semiconductors.co.uk/propiivi5410.htm (all Zn taken from here, 300 K) and above comments are
#define AZnSeWZ 3.98      // c = 6.53, c/a = 1.641,   Zn-Se bond length = 2.45 A, Eg = 
#define AZnSeZB 5.67      //                          Zn-Se bond length = 2.46 A, Eg = 2.8215 eV
#define AZnSWZ  3.811     // c = 6.234, c/a = 1.636,  Zn-S  bond length = 2.34 A, Eg = 3.911  eV
#define AZnSZB  5.421     //                          Zn-S  bond length = 2.34 A, Eg = 3.68   eV
#define AZnTeWZ 4.27      // c = 6.99, c/a = 1.637,   Zn-Te bond length = 2.62 A, Eg = 
#define AZnTeZB 6.10      //                          Zn-Te bond length = 2.64 A, Eg = 2.394  eV 
#define AZnOWZ  3.2495    // c = 5.2069, c/a = 1.602, Zn-O  bond length = 1.95 A, Eg = 3.4    eV

/****************************************************************************/
/* Structure declarations */

typedef struct nanostructure_ {
	int index;
	int nAtoms, nAtomTypes, nMaxBonds;
	int nStackingFaults, nCores;
	char ncType[100], ncCrystalStructure[100];
	char ncAtomSymbol1[4], ncAtomSymbol2[4], centerAtomSymbol[4];
	vector ncSize, ncCenter;
	atom *atoms; // pointer to the first atom
} nanostructure;

/****************************************************************************/
/* Function declarations - public interface */

// Functions that deal with the nanostructure structure - nanostructure.c 


/****************************************************************************/

#endif

/****************************************************************************/
