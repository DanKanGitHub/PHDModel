// Assembles the solution matrix and the RHS

// U_Feet is NNMX3 where each row is the global node number the departure foot is assocaited 
// with and the entires in the row are the x,y coordinates of the departure foot followed by 
// the ellemnt the foot is located in

// A and F come in initialized to zeros

// !!!!!!!!!!!!! addd sparse matrix utilization!!!!!!!!!!!!

#ifndef InitMatZero_H_Guard
#define InitMatZero_H_Guard

#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_IntSerialDenseVector.h"

typedef Epetra_IntSerialDenseVector E_ISDV;
typedef Epetra_IntSerialDenseMatrix E_ISDM;
typedef Epetra_CrsMatrix E_CM;

using namespace std;

void InitMatZero(int Vel_Npe, 
		  int Pre_Npe,
		  int Nem, 
		  int Vel_Nnm,
		  int Pre_Nnm,
		  E_ISDM Vel_Nod, 
		  E_ISDM Pre_Nod,
		  E_ISDM Vel_Nod_BC_Hor,
		  E_ISDM Vel_Nod_BC_Ver,
		  E_ISDM Pre_Nod_BC,
		  E_ISDV All_Proc_Nodes_VP,
		  E_ISDV My_Proc_Eles,
		  int myid,
		  E_CM & A);

#endif