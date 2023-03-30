/****************************************************************************/
/* This file does the printing for the program */

#include "fd.h"

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
