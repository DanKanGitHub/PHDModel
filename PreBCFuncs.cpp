#include "PreBCFuncs.h"
#include "Epetra_SerialDenseVector.h"

typedef Epetra_SerialDenseVector E_SDV;

// essential boundary condition
// Only applies to velocity and all are in global coordinates
double PresEssenEntrance(double x,
			  double y)
{
  
  double z;
  
  z = 12.0;
  return z;
}

// natural boundary conditions
// Pressure
double PreNatLowerPlate(double x,
			double y)
{
  double z;
  
  z = 0.0;
  return z;
}

double PreNatUpperPlate(double x,
		      double y)
{
  double z;
  
  z = 0.0;
  return z;
}

double PreNatEntrance(double x,
		    double y)
{
  double z;
  
  z = 0.0;
  return z;
}

double PreNatExit(double x,
		double y)
{
  double z;
  
  z = 0.0;
  return z;
}