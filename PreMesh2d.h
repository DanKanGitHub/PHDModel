#ifndef PreMesh2d_H_Guard
#define PreMesh2d_H_Guard

#include <cmath>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_IntSerialDenseMatrix E_ISDM;

using std::floor;

void PreMesh2d(double dx, 
		double dy, 
		double XL,
		double YB,
		int NX, 
		int NY,
		int Nem, 
		int & Nnm, 
		int & Npe, 
		E_SDM & Glxy, 
		E_ISDM & Nod,
		E_ISDM & Nod_BC);

#endif