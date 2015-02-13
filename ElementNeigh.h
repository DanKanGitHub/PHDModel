// ElementNeigh determines all the elements which share a node
// The numbering is ddetermined by the local node number, i.e., side 1 is opposite nnode 1.
// This only works for my mesh
// NX is the number of intervals in the horizontal direction.

#ifndef ElementNeigh_H_Guard
#define ElementNeigh_H_Guard

#include "Epetra_IntSerialDenseMatrix.h"

typedef Epetra_IntSerialDenseMatrix E_ISDM;

void ElementNeigh(int Nem,
		  int NX,
		  E_ISDM & Ele_Neigh);

#endif