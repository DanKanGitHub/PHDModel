// Departure Feet algorithm
// Uses only the velocity element and passes back both the departure velocity
// and the departure stress

#ifndef DepartureFoot_H_Guard
#define DepartureFoot_H_Guard

#include "FeetSearch.h"
#include "Shape2d.h"

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_IntSerialDenseMatrix E_ISDM;
typedef Epetra_SerialDenseVector E_SDV;

void DepartureFoot(int Vel_Flag, 
		   int Nem, 
		   int Vel_Nnm,
		   int NX,
		   double TIMESTEP,
		   double dx,
		   double dy,
		   double *Vel, 
		   E_SDM Vel_Glxy, 
		   int Vel_Npe, 
		   E_ISDM Vel_Nod, 
		   E_ISDM Ele_Neigh,
		   E_SDV Ini_Foot,
		   E_SDV & Depart_Foot,
		   int & New_Ele); 

#endif