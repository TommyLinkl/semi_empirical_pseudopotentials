/*****************************************************************************/
// This file handles errors: primarily just writes error messages and exits
// the programs

#include "fd.h"

/*****************************************************************************/

void nerror(char *str) {
  FILE *pf;

  pf = fopen("error", "w");
  fprintf(pf, "%s\n", str);
  
  exit(EXIT_FAILURE);
}

/*****************************************************************************/
// Write error message to the file pointed to by pf (often will be an error.dat,
// stderr or stdout) and exit the program 

void memoryError(char *str, FILE *pf) {

	writeSeparation(pf);
	fprintf(pf, "MEMORY ERROR:\n\n");
	fprintf(pf, "Could not allocate memory for %s\n\n", str);
	fprintf(pf, "The program is exiting due to this memory allocation failure\n");
	writeSeparation(pf);
	
	exit(EXIT_FAILURE);
}

/*****************************************************************************/
