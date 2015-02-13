// The initial velocity profile has zero vertical component and quadratic horizontal component with zero 
// on the boundaries.

#ifndef InitialBio_H_Guard
#define InitialBio_H_Guard

#include "Epetra_SerialDenseMatrix.h"

#include <cmath>

typedef Epetra_SerialDenseMatrix E_SDM;

void InitialBio(int NX,
		E_SDM Glxy,
		double *Bio);

#endif