#ifndef StrNodeDepartureFeet_H_Gaurd
#define StrNodeDepartureFeet_H_Gaurd

#include "DepartureFoot.h"

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_IntSerialDenseMatrix E_ISDM;
typedef Epetra_SerialDenseVector E_SDV;
typedef Epetra_IntSerialDenseVector E_ISDV;

void StrNodeDepartureFeet(int myid,
			  int Nem,
			  int Vel_Npe,
			  int Str_Npe,
			  int VEL_FLAG,
			  int Vel_Nnm,
			  int NX,
			  int STRESS_FLAG,
			  double TIMESTEP,
			  double dx,
			  double dy,
			  E_SDM Vel_Glxy,
			  E_ISDM Vel_Nod,
			  E_ISDM Ele_Neigh,
			  E_ISDV All_Proc_Nodes,
			  E_ISDV My_Proc_Eles,
			  double *Vel_Old,
			  E_SDM & StrNodeDepartFootx,
			  E_SDM & StrNodeDepartFooty,
			  E_ISDM & StrNodeDepartElement,
		  E_ISDM Str_Nod,
		  E_SDM Str_Glxy,
		  E_SDM Str, 
		  E_SDM Str_Old);
#endif