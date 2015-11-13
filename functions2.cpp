#include "functions2.h"

// essential boundary condition
// Only applies to velocity and all are in global coordinates
double PresEssenEntrance(double x,
			  double y)
{
  
  double z;
  
  z = 0.0; // 12.0;
  return z;
}

void EssenEntrance(double x,
		   double y,
		   double Cur_Time,
		   double Timestep,
		   E_SDV & z)
{
  double ramp, numRampDownTimeStep;
  double freq, numTimeStepRampChange, rampTime, numRampUpTimeStep;
  double pi = 3.14159265358979323846;
  
  // The number of time steps that the inflow will take to change from one condition to the next.
  // Set this so that it is equal to 5 to 10 time steps.
  numTimeStepRampChange = 5.0;
  
  // The time it will take for the change in inflow to take place
  rampTime = numTimeStepRampChange * Timestep;
  
  // The effect on the period of the trig functions used for ramping
  // up and down.
  freq = pi / (2.0 * rampTime);
  
  // Number of time steps before the inflow is turned off
  numRampDownTimeStep = 200.0; // approx 10 seconds
  
  // Number of time steps before the inflow is turned on
  // Not used at this time.
  numRampUpTimeStep = 0.0;

  // Ramping up
  if(Cur_Time < rampTime) {
    ramp = sin(freq * Cur_Time); // ramp is the rate at which the quadratic boundary function increases to 1
    z(0) = -4 * ramp * y * (y - 1.0);
    
    // Running at constant inflow
  } else if(((Cur_Time >= rampTime) && (Cur_Time <= numRampDownTimeStep * Timestep))) {
    z(0) = -4.0 * y * (y - 1.0); // Coeff if 
    
    // Ramping down
  } else if((Cur_Time > numRampDownTimeStep * Timestep) && (Cur_Time < (numRampDownTimeStep + numTimeStepRampChange) * Timestep)) { // Approximately 6 seconds
    ramp = cos(freq * (Cur_Time - numRampDownTimeStep * Timestep));
    z(0) = -4 * ramp * y * (y - 1.0);
    
    // Inflow turned off
  } else {
    z(0) = 0.0;
  };
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

void BioEssenExit(double x,
	       double y,
	       double *Bio,
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