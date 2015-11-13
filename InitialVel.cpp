// The initial velocity profile has zero vertical component and quadratic horizontal component with zero 
// on the boundaries.

#include "InitialVel.h"

void InitialVel(int Nnm,
		E_SDM Glxy,
		double *Vel) 
{
  double x, y;
  
  for (int i = 0; i <= Nnm-1; i++)
  {
    x = Glxy(i,0);
    y = Glxy(i,1);
    Vel[i] = 0.0; //HorizVel(x,y);
  }
}