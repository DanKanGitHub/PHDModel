// Assembles the solution matrix and the RHS

// U_Feet is NNMX3 where each row is the global node number the departure foot is assocaited 
// with and the entires in the row are the x,y coordinates of the departure foot followed by 
// the ellemnt the foot is located in

// A and F come in initialized to zeros

// !!!!!!!!!!!!! addd sparse matrix utilization!!!!!!!!!!!!

#ifndef InitMatZeroBio_H_Guard
#define InitMatZeroBio_H_Guard

#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_IntSerialDenseVector.h"

typedef Epetra_IntSerialDenseVector E_ISDV;
typedef Epetra_IntSerialDenseMatrix E_ISDM;
typedef Epetra_CrsMatrix E_CM;

using namespace std;

void InitMatZeroBio(int Vel_Npe, 
		    int Nem,
		    int Vel_Nnm,
		    E_ISDM Vel_Nod,
		    E_ISDM Bio_Nod_BC,
		    E_ISDV All_Proc_Nodes_Bio,
		    E_ISDV My_Proc_Eles,
		    int myid,
		    E_CM & A);

#endif