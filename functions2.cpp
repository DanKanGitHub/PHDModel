#include "functions2.h"
// #include "Epetra_SerialDenseVector.h"

// typedef Epetra_SerialDenseVector E_SDV;

// essential boundary condition
// Only applies to velocity and all are in global coordinates
double PresEssenEntrance(double x,
			  double y)
{
  
  double z;
  
  z = 12.0;
  return z;
}

void EssenEntrance(double x,
		   double y,
		   E_SDV & z)
{
  z(0) = -4.0 * y * (y - 1.0); // Coeff if 
  z(1) = 0.0;
}

void EssenUpperPlate(double x,
		     double y,
		     E_SDV & z)
{
  z(0) = 0.0;
  z(1) = 0.0;
}

void EssenLowerPlate(double x,
		     double y,
		     E_SDV & z)
{
  z(0) = 0.0;
  z(1) = 0.0;
}

void EssenExit(double x,
	       double y,
	       E_SDV & z)
{
  z(0) = 0.0;
  z(1) = 0.0;
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

// Velocity
void VelNatLowerPlate(double x,
		      double y,
		      E_SDV & z)
{
  z(0) = 0.0;
  z(1) = 0.0;
}

void VelNatUpperPlate(double x,
		      double y,
		      E_SDV & z)
{
  z(0) = 0.0;
  z(1) = 0.0;
}

void VelNatEntrance(double x,
		    double y,
		    E_SDV & z)
{
  z(0) = 0.0;
  z(1) = 0.0;
}

void VelNatExit(double x,
		double y,
		E_SDV & z)
{
  z(0) = 0.0;
  z(1) = 0.0;
}

// function return 1-D shape function at natural coordinate xi, order 1 or
// 2, page 346
void Sf1D(double xi,
	  int order,
	  E_SDV & z)
{
  if (order == 1)
  {
    z(0) = 0.5 * (1 - xi);
    z(1) = 0.5 * (1 + xi);
  }
  else if (order == 2)
  {
    z(0) = -0.5 * xi * (1 - xi);
    z(1) = 1.0 - xi * xi;
    z(2) = 0.5 * xi * (1 + xi);
  }
}