// conformation tensor solver.  Uses linear elements only
// Can reduce the number of Gauss points for this program and have the same accuracy

#ifndef Conformation_H_Guard
#define Conformation_H_Guard

// #include <cmath>
// #include "Shape2d.h"
// #include "BioWeightFunc.h"
// #include "RetardationDividedByRelaxation.h"
#include "NodeStress.h"
#include "DepartureFoot.h"

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_IntSerialDenseMatrix E_ISDM;
typedef Epetra_SerialDenseVector E_SDV;
typedef Epetra_IntSerialDenseVector E_ISDV;

void Conformation(double TIMESTEP, 
		  double Sol_Vis,
		  double Poly_vis,
		  int NX,
		  double RETARD_TIME,
		  double Sol_Density,
		  double Poly_Density,
		  double T_ZERO,
		  double L_ZERO,
		  int Vel_Npe, 
		  int Str_Npe, 
		  int Nem, 
		  int Vel_Nnm,
		  int myid,
		  E_ISDM Vel_Nod, 
		  E_ISDM Str_Nod,
		  E_SDM Vel_Glxy, 
		  E_SDM Str_Glxy,
		  E_ISDM Ele_Neigh,
		  int VEL_FLAG, 
		  int STRESS_FLAG, 
		  double *Vel, 
		  double *Vel_Old,
		  double *Bio,
		  E_SDM Str, 
		  E_SDM Str_Old,
		  E_SDM StrNodeDepartFootx,
		  E_SDM StrNodeDepartFooty,
		  E_ISDM StrNodeDepartElement,
		  E_ISDV All_Proc_Nodes,
		  E_ISDV My_Proc_Eles,
		  E_SDM & Str_New);

#endif