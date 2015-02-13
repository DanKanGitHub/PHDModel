// The initial velocity profile has zero vertical component and quadratic horizontal component with zero 
// on the boundaries.

#include "InitialVel.h"

void InitialVel(int Nnm,
		E_SDM Glxy,
		double *Vel) 
{
  double y;
  for (int i = 0; i <= Nnm-1; i++)
  {
    y = Glxy(i,1);
    Vel[i] = -3.0 * y * (y - 1.0);
//     Vel(i + Nnm) = 0.0;
  }
}