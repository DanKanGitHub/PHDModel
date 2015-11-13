// Velocity and its partials called VelFunctions
// Used by InitialVel.cpp

#include "VelFunctions.h"

double HorizVel(double x,
		double y)
{
  double z;
  
  z = 0.0; // -4.0 * y * (y - 1.0);// The max vel is 1 with a coeff of -4.
  
  return z;
}

double VertVel(double x,
	      double y)
{
  double z;
  
  z = 0.0;
  
  return z;
}

double PartialHorizVelPartialx(double x,
			      double y)
{
  double z;
  
  z = 0.0;
  
  return z;
}

double PartialHorizVelPartialy(double x,
			      double y)
{
  double z;
  
  z = 0.0; // -8.0 * y + 4.0;
  
  return z;
}

double PartialVertVelPartialx(double x,
			    double y)
{
  double z;
  
  z = 0.0;
  
  return z;
}

double PartialVertVelPartialy(double x,
			    double y)
{
  double z;
  
  z = 0.0;
  
  return z;
}