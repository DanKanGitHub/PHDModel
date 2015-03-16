// Viscosity * 0.5 * (Gradu + GRaduTrans)
// The other values not discribed here are initialized to zero.

#ifndef InitialStr_H_Guard
#define InitialStr_H_Guard

#include "BioWeightFunc.h"
#include "RetardationDividedByRelaxation.h"
#include "VelFunctions.h"

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_IntSerialDenseMatrix E_ISDM;

void InitialStr(int Nem, 
		int Str_Npe,
		double Sol_Vis,
		double Poly_Vis,
		double Sol_Density,
		double Poly_Density,
		double L_ZERO,
		double T_ZERO,
		E_ISDM Vel_Nod,
		E_ISDM Str_Nod,
		E_SDM Glxy, 
		double *Bio,
		E_SDM & Str);

#endif