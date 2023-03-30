/****************************************************************************/
//
//
//
/****************************************************************************/

#ifndef EH_H
#define EH_H

/****************************************************************************/
/* These are the library functions that are used */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>

/****************************************************************************/
/* Macro definitions: common unit conversion and multiplication schemes
  that help make the program more readable */



/****************************************************************************/
/* Structure declarations */



/****************************************************************************/
/* Function declarations - public interface */

/* Functions that handle errors/ perform error checks - errorHandling.c */
void memoryError(const char *errorMessage);

/****************************************************************************/

#endif

/****************************************************************************/
