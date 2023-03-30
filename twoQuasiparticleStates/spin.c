/*****************************************************************************/
//
// This file deals with the spin quantum numbers of multiquasiparticle states 
//
/*****************************************************************************/

#include "qp.h"

/*****************************************************************************/
// Returns the spin z projection quantum number, Sz, for a noninteracting two qp state
// Sz = Sz1 + Sz2

double retSpinZProjNonintTwoQPState(nonintTwoQPState nonintTwoQP) {
	double Sz;

	// Sz = Sz1 + Sz2
	Sz = (nonintTwoQP.qp1->spinZ + nonintTwoQP.qp2->spinZ);

	return Sz;
}

/*****************************************************************************/
// Returns the total spin quantum number, S, for a noninteracting two qp state
// <S^2> = <Psi1Psi2|S^2|Psi1Psi2> is proportional to S(S+1)
// S^2 = S1^2 + S2^2 + 2*S12*S12 = S1^2 + S2^2 + S1+S2- + S1-S2+ + 2S1zS2z
// <S^2> = <Psi1|S1|Psi1> + <Psi2|S2|Psi2> + 0 + 0 + 2*S1z*S2z
// <S^2> = <Psi1|S1|Psi1> + <Psi2|S2|Psi2> + 2*S1z*S2z
// <S^2> = S1*(S1+1) + S2*(S2+1) + 2*S1z*S2z

double retTotalSpinNonintTwoQPState(nonintTwoQPState nonintTwoQP) {
	double S;

	// <S^2> = <S1^2> + <S2^2> + 2*<S1z>*<S2z> 
	// <S^2> = S1*(S1+1) + S2*(S2+1) + 2*S1z*S2z  
	S  = ( nonintTwoQP.qp1->spin * (nonintTwoQP.qp1->spin + 1.0) );
	S += ( nonintTwoQP.qp2->spin * (nonintTwoQP.qp2->spin + 1.0) );
	S += ( 2.0 * nonintTwoQP.qp1->spinZ * nonintTwoQP.qp2->spinZ ); 

	if ( fabs(S - 2.0) < EPS) {
		return 1.0; // s = 1 for <S^2> == 2
	}
	else {
		return 0.0; // s = 0 for <S^2> == 0
	}

}

/*****************************************************************************/
// Returns the z projected spin quantum number, Sz, for an interacting two qp state
// Sz = SumOverRS(Crs*Szr*Szs)

double retSpinZProjIntTwoQPState(intTwoQPState intTwoQP, long nNonintTwoQPStates) {
	long i;
	double Sz = 0.0;
	nonintTwoQPState *nonintTwoQP = intTwoQP.niTwoQP;

	// Sz = SumOverRS(Crs*Szr*Szs)
	for (i = 0; i < nNonintTwoQPStates; i++) {
		Sz += ( intTwoQP.Crs[i] * ( nonintTwoQP[i].qp1->spinZ + nonintTwoQP[i].qp2->spinZ) ); 
	}

	return Sz;
}

/*****************************************************************************/
// Returns the total spin quantum number, S, for an interacting two qp state
// <S^2> = SumOverRS(Crs*<Srs^2>)

double retTotalSpinIntTwoQPState(intTwoQPState intTwoQP, long nNonintTwoQPStates) {
	long i;
	double S = 0.0;
	nonintTwoQPState *nonintTwoQP = intTwoQP.niTwoQP;

	// <S^2> = SumOverRS(Crs*<Srs^2>)
	for (i = 0; i < nNonintTwoQPStates; i++) {
		S += ( intTwoQP.Crs[i] * retTotalSpinNonintTwoQPState(nonintTwoQP[i]) ); 
	}

	return S;
}

/*****************************************************************************/

