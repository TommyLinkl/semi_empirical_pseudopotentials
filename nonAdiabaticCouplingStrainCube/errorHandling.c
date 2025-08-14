/*****************************************************************************/

#include "qp.h"

/*****************************************************************************/
void memoryError(char *str) {
  FILE *pf;

  pf = fopen("error.dat", "w");
  fprintf(pf, "%s\n", str);
  fclose(pf);

  writeSeparation(stdout);
  writeCurrentTime(stdout);
  fprintf(stdout, "There was a memory allocation error that caused the program to exit!\n");
  fprintf(stdout, "Error message: %s\n", str);
  writeSeparation(stdout);
  fflush(stdout);
  
  exit(EXIT_FAILURE);
}

/*****************************************************************************/
void nerror(char *str)
{
  FILE *pf;

  pf = fopen("error","w");
  fprintf(pf,"%s\n",str);
  exit(1);
}

/*****************************************************************************/