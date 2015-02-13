// The initial velocity profile has zero vertical component and quadratic horizontal component with zero 
// on the boundaries.

#include "InitialPre.h"

void InitialPre(int Vel_Nnm,
		double *Pre) 
{
  Pre[2 * Vel_Nnm] = 0.0;
}