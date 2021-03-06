// Assembles the solution matrix and the RHS

// U_Feet is NNMX3 where each row is the global node number the departure foot is assocaited 
// with and the entires in the row are the x,y coordinates of the departure foot followed by 
// the ellemnt the foot is located in

// A and F come in initialized to zeros

// !!!!!!!!!!!!! addd sparse matrix utilization!!!!!!!!!!!!

#ifndef EquibSparseAssembly_H_Guard
#define EquibSparseAssembly_H_Guard

#include "Shape2d.h"
#include "VelEssenBoundary2d.h"
#include "PreEssenBoundary2d.h"
#include "NatBoundary2d.h"
#include "BioWeightFunc.h"
#include "RetardationDividedByRelaxation.h"

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

void EquibSparseAssembly(double Sol_Vis, 
			  double Poly_Vis, 
			  double Sol_Density,
			  double Poly_Density,
			  double T_ZERO,
			  double L_ZERO,
			  double TIMESTEP, 
			  int NX,
			  int NY,
			  int NGP,
			  int Vel_Npe, 
			  int Pre_Npe,
			  int Nem, 
			  int N_TRI_QUAD,
			  int Vel_Nnm,
			  int Pre_Nnm,
			  E_ISDM Vel_Nod,
			  E_ISDM Pre_Nod,
			  E_ISDM Vel_Nod_BC_Hor,
			  E_ISDM Vel_Nod_BC_Ver,
			  E_ISDM Pre_Nod_BC, 
			  E_SDM Vel_Glxy, 
			  E_SDM Pre_Glxy, 
			  E_ISDM Ele_Neigh,
			  int VEL_FLAG,
			  int PRE_FLAG, 
			  E_SDM Tri_Quad_Pt, 
			  E_SDV Tri_Quad_Wt,
			  E_SDM GAUSPT,
			  E_SDM GAUSWT,
			  E_SDM Str, 
			  double *Bio,
			  double *Vel_Old,
			  E_ISDV All_Proc_Nodes_VP,
			  E_ISDV My_Proc_Eles,
      // 		    E_SDV Shared_Nodes_VP,
			  int myid,
			  E_CM & A, 
			  E_V & b);
#endif