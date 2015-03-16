// The initial velocity profile has zero vertical component and quadratic horizontal component with zero 
// on the boundaries.

#include "InitialVel.h"

void InitialVel(int Nnm,
		E_SDM Glxy,
		double *Vel) 
{
  double x, y;
  
  x = 0.0;
  
  for (int i = 0; i <= Nnm-1; i++)
  {
    y = Glxy(i,1);
    Vel[i] = HorizVel(x,
		      y);
  }
}