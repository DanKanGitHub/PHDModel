// function to impose boundary conditions (Essential and Natural) for 2D finite element code

#ifndef VelBoundary2d_H_Guard
#define VelBoundary2d_H_Guard

#include "functions2.h"
// #include "Epetra_SerialDenseMatrix.h"
// #include "Epetra_SerialDenseVector.h"

// typedef Epetra_SerialDenseMatrix E_SDM;
// typedef Epetra_SerialDenseVector E_SDV;

void VelEssenBoundary2d(int NX,
			int Vel_Nnm,
			int ii, 
			double x,
			double y,
			E_SDV & Vel_Gn);
#endif
