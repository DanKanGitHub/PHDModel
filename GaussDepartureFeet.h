#ifndef GaussDepartureFeet_H_Gaurd
#define GaussDepartureFeet_H_Gaurd

#include "DepartureFoot.h"

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_IntSerialDenseMatrix E_ISDM;
typedef Epetra_SerialDenseVector E_SDV;
typedef Epetra_IntSerialDenseVector E_ISDV;

void GaussDepartureFeet(int myid,
			int Nem,
			int N_TRI_QUAD,
			int Vel_Npe,
			int VEL_FLAG,
			int Vel_Nnm,
			double TIMESTEP,
			E_SDM Vel_Glxy,
			E_SDM Tri_Quad_Pt,
			E_ISDM Vel_Nod,
			E_ISDM Ele_Neigh,
			E_ISDV My_Proc_Eles,
			double *Vel_Old,
			E_SDM & GaussDepartFootx,
			E_SDM & GaussDepartFooty,
			E_ISDM & DepartureElement);
#endif