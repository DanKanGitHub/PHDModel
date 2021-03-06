// Assembles the solution matrix and the RHS

// U_Feet is NNMX3 where each row is the global node number the departure foot is assocaited 
// with and the entires in the row are the x,y coordinates of the departure foot followed by 
// the ellemnt the foot is located in

// A and F come in initialized to zeros

// !!!!!!!!!!!!! addd sparse matrix utilization!!!!!!!!!!!!

#ifndef BioAssembly_H_Guard
#define BioAssembly_H_Guard

#include "Shape2d.h"

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_IntSerialDenseMatrix E_ISDM;
typedef Epetra_SerialDenseVector E_SDV;
typedef Epetra_IntSerialDenseVector E_ISDV;

typedef Epetra_CrsMatrix E_CM;
typedef Epetra_Vector E_V;

using namespace std;

void BioSparseAssembly(double DIFF_COEFF,
			double TIMESTEP,
			int NX,
			int NY,
			int NGP,
			int Vel_Npe,
			int Nem, 
			int N_TRI_QUAD,
			int Vel_Nnm,
			E_ISDM Vel_Nod, 
			E_ISDM Bio_Nod_BC,
			E_SDM Vel_Glxy,
			E_SDM Tri_Quad_Pt, 
			E_SDV Tri_Quad_Wt,
			E_SDM GAUSPT,
			E_SDM GAUSWT,
			double *Vel,
			double *Vel_Old,
			double *Bio_Old,
			E_ISDV All_Proc_Nodes_Bio,
			E_ISDV My_Proc_Eles,
// 			E_SDV Shared_Nodes_Bio,
			int myid,
			E_CM & A_Bio, 
			E_V & b_Bio);
#endif