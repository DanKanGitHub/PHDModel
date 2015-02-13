// Determines which element the point x is in.
// CEle:	current element number
// NEIGH:	All elements which are the Neighbor of CEle, size 4XNEM

#ifndef FeetSearch_H_Guard
#define FeetSearch_H_Guard

#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

typedef Epetra_SerialDenseVector E_SDV;
typedef Epetra_IntSerialDenseMatrix E_ISDM;

int FeetSearch(E_SDV x, 
	       int Cur_Ele, 
	       E_ISDM Ele_Neigh);

#endif