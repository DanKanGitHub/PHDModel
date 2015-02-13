// function to impose boundary conditions (Essential and Natural) for 2D finite element code

#include "PreEssenBoundary2d.h"

double PreEssenBoundary2d(double x,
			  double y)
{
  // Variable declarations
  double Pre_Gn;

  // Essential boundary conditions
  Pre_Gn = PresEssenEntrance(x,y); // could be set to 12.0
  
  return Pre_Gn;
}