// Stress mesh2d builds linear for Str_Flag = 1 and quadratic for Str_Flag = 2
// Linear follows pressure and Quadratic follows velocity and always w/o BCs

#include "StrMesh2d.h"

void StrMesh2d(double dx, 
		  double dy, 
		  double XL,
		  double YB,
		  int Str_Flag,
		  int NX, 
		  int NY,
		  int Nem, 
		  int & Nnm, 
		  int & Npe, 
		  E_SDM & Glxy, 
		  E_ISDM & Nod)
{
  // Variable declarations
  int Eln1, Eln2;
  
  if(Str_Flag == 1) // Linear Case
  {
    // % size of various structures
    Nnm = (NY + 1) * (NX + 1);
    Npe = 3;
    
    Glxy.Shape(Nnm,2);
    for(int i = 0; i <= Nnm - 1; i++)
    {
      Glxy(i,0) = 0.0;
      Glxy(i,1) = 0.0;
    }
    
    Nod.Shape(Nem,Npe); // Added a dimension for the EBC, NBC flag.

    // global coordinate of nodes
    for (int i = 0; i <= Nnm-1; i++)
    {
      Glxy(i,0) = (i % (NX + 1)) * dx + XL;
      Glxy(i,1) = floor(i / (NX + 1)) * dy + YB;
    }

    // connectivity matrix
    for (int j = 0; j <= NY-1; j++)	// indexes the horizontal elements
    {
      for (int i = 0; i <= NX-1; i++)	// indexes the vertical elements
      {
	Eln1 = j*2*NX + 2*(i + 1) - 1;	// lower element number
	Eln2 = Eln1 + 1;		// upper element number

	Nod(Eln1-1,0) = j*(NX+1) + i + 1;
	Nod(Eln1-1,1) = Nod(Eln1-1,0) + 1;
	Nod(Eln1-1,2) = Nod(Eln1-1,1) + NX + 1;

	Nod(Eln2-1,0) = Nod(Eln1-1,0);
	Nod(Eln2-1,1) = Nod(Eln1-1,2);
	Nod(Eln2-1,2) = Nod(Eln2-1,1) - 1;
      }
    }
  }
  else // Quadratic case
  {
    
    NX = 2 * NX;
    NY = 2 * NY;
    
    dx = dx / 2;
    dy = dy / 2;
    
    // size of various structures
    Nnm = (NY + 1) * (NX + 1);
    Npe = 6;
    
    Glxy.Shape(Nnm,2);
    for(int i = 0; i <= Nnm - 1; i++)
    {
      Glxy(i,0) = 0.0;
      Glxy(i,1) = 0.0;
    }
    
    Nod.Shape(Nem,Npe);

    // global coordinate of nodes
    for (int i = 0; i <= Nnm-1; i++)
    {
      Glxy(i,0) = (i % (NX + 1)) * dx + XL;
      Glxy(i,1) = floor(i / (NX + 1)) * dy + YB;
    }

    // connectivity matrix
    for (int j = 0; j <= NY/2-1; j++) 	// indexes the horizontal elements
    {
      for (int i = 0; i <= NX/2-1; i++)// indexes the vertical elements
      {
	Eln1 = j*NX + 2*(i + 1) - 1;	// lower element number,
	Eln2 = Eln1 + 1;		// upper element number

	Nod(Eln1-1,0) = 2 * j * (NX + 1) + 2 * i + 1; //Global node number
	Nod(Eln1-1,1) = Nod(Eln1-1,0) + 2;
	Nod(Eln1-1,2) = Nod(Eln1-1,1) + 2*(NX + 1);
	Nod(Eln1-1,3) = Nod(Eln1-1,0) + 1;
	Nod(Eln1-1,4) = Nod(Eln1-1,1) + NX + 1;
	Nod(Eln1-1,5) = Nod(Eln1-1,4) - 1;

	Nod(Eln2-1,0) = Nod(Eln1-1,0);
	Nod(Eln2-1,1) = Nod(Eln1-1,2);
	Nod(Eln2-1,2) = Nod(Eln2-1,1) - 2;
	Nod(Eln2-1,3) = Nod(Eln1-1,5);
	Nod(Eln2-1,4) = Nod(Eln2-1,1) - 1;
	Nod(Eln2-1,5) = Nod(Eln2-1,3) - 1;

      }
    } 
  }
}