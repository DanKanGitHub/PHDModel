#ifndef InitialGuess_H_Guard
#define InitialGuess_H_Guard

#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_Vector.h"

typedef Epetra_IntSerialDenseMatrix E_ISDM;
typedef Epetra_Vector E_V;

void InitialGuess(int *ProcNodes, 
		  int Npe,
		  int Vel_Nnm,
		  int *ProcEles,
		  E_ISDM Vel_Nod,
		  E_ISDM Pre_Nod,
		  int TotalNumEles,
		  int myid,
		  double *PreviousSoln, 
		  E_V & x_VP);

#endif