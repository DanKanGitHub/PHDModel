#ifndef StrMesh2d_H_
#define StrMesh2d_H_

#include <cmath>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_IntSerialDenseMatrix E_ISDM;

using std::floor;

void StrMesh2d(double dx, 
		  double dy, 
		  double XL,
		  double YB,
		  int Str_Flag,
		  int NX, 
		  int NY,
		  int Nem, 
		  int & Nnm, 
		  int & Npe, 
		  E_SDM & Glxy, 
		  E_ISDM & Nod);

#endif