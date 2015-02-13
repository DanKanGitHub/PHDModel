// The initial velocity profile has zero vertical component and quadratic horizontal component with zero 
// on the boundaries.

#ifndef InitialVel_H_Guard
#define InitialVel_H_Guard

#include "Epetra_SerialDenseMatrix.h"

typedef Epetra_SerialDenseMatrix E_SDM;

void InitialVel(int Nnm,
		E_SDM Glxy,
		double *Vel);

#endif