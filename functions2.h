// essential boundary condition

#ifndef FUNCTIONS2_H_Guard
#define FUNCTIONS2_H_Guard

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

typedef Epetra_SerialDenseVector E_SDV;
typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_IntSerialDenseMatrix E_ISDM;

// pressure Essential boundary conditions
double PresEssenEntrance(double x,
			  double y);

// Velocity Essential boundary conditions
void EssenEntrance(double x,
		   double y,
		   E_SDV & z);

void EssenUpperPlate(double x,
		     double y,
		     E_SDV & z);

void EssenLowerPlate(double x,
		     double y,
		     E_SDV & z);

void EssenExit(double x,
	       double y,
	       E_SDV & z);

// natural boundary conditions
// Pressure
double PreNatLowerPlate(double x,
			double y);

double PreNatUpperPlate(double x,
		      double y);

double PreNatEntrance(double x,
		    double y);

double PreNatExit(double x,
		double y);

// Velocity
void VelNatLowerPlate(double x,
		      double y,
		      E_SDV & z);

void VelNatUpperPlate(double x,
		      double y,
		      E_SDV & z);

void VelNatEntrance(double x,
		    double y,
		    E_SDV & z);

void VelNatExit(double x,
		double y,
		E_SDV & z);

void Sf1D(double xi,
	  int order,
	  E_SDV & z);

#endif