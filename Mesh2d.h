#ifndef Mesh2d_H_Guard
#define Mesh2d_H_Guard

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_IntSerialDenseMatrix E_ISDM;

void Mesh2d(double dx, 
	    double dy, 
	    double XL, 
	    double XR, 
	    double YB, 
	    double YT, 
	    int NX, 
	    int NY, 
	    int FLAG,
	    int Nem, 
	    int & Nnm, 
	    int & Npe, 
	    E_SDM & Glxy, 
	    E_ISDM & Nod,
	    E_ISDM & Nod_BC_Hor,
	    E_ISDM & Nod_BC_Ver);

#endif