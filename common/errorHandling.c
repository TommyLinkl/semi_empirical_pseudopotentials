/****************************************************************************/
//
// This file contains functions related to handling errors
//
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nc.h"

/****************************************************************************/
// Prints error message and exits the program

void memoryError(const char *errorMessage) {
	
	// Write error message
	writeSeparation(stdout);
	fprintf(stdout, "MEMORY ALLOCATION ERROR:\n");
	writeCurrentTime(stdout);
	fprintf(stdout, "\n%s\n", errorMessage);
	fprintf(stdout, "The program is exiting due to this memory allocation failure\n");
	writeSeparation(stdout);
	fflush(stdout);

	// Exit the program
	exit(EXIT_FAILURE);
}

/****************************************************************************/