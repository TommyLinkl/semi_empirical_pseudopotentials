/****************************************************************************/
/* This file does the printing for the program */

#include "fd.h"

/****************************************************************************/
// Writes a separation then the current time followed by the message 
// Used to simplify timing writing in the begging of a function
 
void writeFunctionStartTime(char *message, FILE *pf) {
	writeSeparation(pf);
	writeCurrentTime(pf);
	fprintf(pf, "\n%s\n", message);

	return;
}

/****************************************************************************/
// prints out current time to stdout 
 
void writeCurrentTime(FILE *pf) {
	time_t startTime;

	startTime = time(NULL);
	fprintf(pf, ctime(&startTime));

	return;
}

/****************************************************************************/

void writeSeparation(FILE *pf) {
	fprintf(pf, "\n******************************************************************************\n\n");
  
	return;
}

/****************************************************************************/
