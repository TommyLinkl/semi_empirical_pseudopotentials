/****************************************************************************/
//
//
//
/****************************************************************************/

#include "fd.h"

/****************************************************************************/
// Returns the number of converged eigenstates with energies between eMin and
// eMax that have sige = sqrt(<E^2>-<E>^2) less than DEPS.
// DEPS = 1e-3 in fd.h

long num_eigenstates_energy_range(double *eval, double *sige, double eMin, double eMax, long maxNumStates) {
	long i, num_eigenstates = 0;

	for (i = 0; i < maxNumStates; i++) {
		if (eval[i] > eMin && eval[i] < eMax && sige[i] < DEPS) num_eigenstates++;
	}

	return (num_eigenstates);
}

/****************************************************************************/
// Stores the indexes of the eigenstates with energies between eMin and eMax
// with sige = sqrt(<E^2>-<E>^2) less than DEPS in eigenstate_index_list 
// and returns the number of eigenstates found 

long get_eigenstate_list(long *eigenstate_index_list, double *eval, double *sige, double eMin, double eMax, long maxNumStates) {
	long i, count = 0; 

	for (i = 0; i < maxNumStates; i++) {
		if (eval[i] > eMin && eval[i] < eMax && sige[i] < DEPS) {
			eigenstate_index_list[count] = i;
			count++;
		}
	}

	return (count);
}

/****************************************************************************/
