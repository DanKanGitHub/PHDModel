// essential boundary condition

#ifndef PREBCFUNCS_H_Guard
#define PREBCFUNCS_H_Guard

// pressure Essential boundary conditions
double PresEssenEntrance(double x,
			  double y);

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

#endif