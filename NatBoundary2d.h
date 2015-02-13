// function to impose boundary conditions (Essential and Natural) for 2D finite element code

#ifndef NatBoundary2d_H_Guard
#define NatBoundary2d_H_Guard

#include "functions2.h"
// #include "Epetra_SerialDenseMatrix.h"
// #include "Epetra_IntSerialDenseMatrix.h"
// #include "Epetra_SerialDenseVector.h"
// 
// typedef Epetra_SerialDenseMatrix E_SDM;
// typedef Epetra_IntSerialDenseMatrix E_ISDM;
// typedef Epetra_SerialDenseVector E_SDV;

void NatBoundary2d(E_SDM GAUSPT, 
		    E_SDM GAUSWT, 
		    int NX,
		    int Ngp,
		    E_SDM Glxy, 
		    E_ISDM Nod,
		    int Gbe,
		    E_SDV & Vel_Gn);
#endif
