// function to impose boundary conditions (Essential and Natural) for 2D finite element code

//  FC and AC enter initialized to all zeros and to the correct size.

#include "VelEssenBoundary2d.h"

void VelEssenBoundary2d(int NX,
			int Vel_Nnm,
			int ii, 
			double x,
			double y,
			E_SDV & Vel_Gn)
{
  // Essential boundary conditions
  // If on the lower plate
  // ii+1 to compensate for the - 1 above
  if(ii + 1 <= 2 * NX + 1)
  {
    EssenLowerPlate(x,y,Vel_Gn); // Note recycling the Vel_Gn variable
  }
// If on the upper plate
  else if(ii+1 > Vel_Nnm - (2 * NX + 1))
  {
    EssenUpperPlate(x,y,Vel_Gn);
  }
  else
  {
    EssenEntrance(x,y,Vel_Gn);
  }
}