// function to impose boundary conditions (Essential and Natural) for 2D finite element code

//  FC and AC enter initialized to all zeros and to the correct size.

#include "NatBoundary2d.h"
// #include "Epetra_SerialDenseMatrix.h"
// #include "Epetra_IntSerialDenseMatrix.h"
// #include "Epetra_SerialDenseVector.h"

// typedef Epetra_SerialDenseMatrix E_SDM;
// typedef Epetra_IntSerialDenseMatrix E_ISDM;
// typedef Epetra_SerialDenseVector E_SDV;

void NatBoundary2d(E_SDM GAUSPT, 
		    E_SDM GAUSWT, 
		    int NX,
		    int Ngp,
		    E_SDM Glxy, 
		    E_ISDM Nod,
		    int Gbe,
		    E_SDV & Vel_Gn)
{
  // Variable declarations
  E_SDM Lin_Tri_En, Quad_Tri_En;
  E_SDV Sf;
  int Ln1, Ln2, jj, NodLn1, NodLn2, Lbed;
  double Xi, Const, x , y, s, x1, y1, x2, y2, h, temp1, temp2;
  
  // edge-node connectivity 
  // linear triangle
  Lin_Tri_En.Shape(3,2);
  Lin_Tri_En(0,0) = 1;
  Lin_Tri_En(0,1) = 2;
  Lin_Tri_En(1,0) = 2;
  Lin_Tri_En(1,1) = 3;
  Lin_Tri_En(2,0) = 3;
  Lin_Tri_En(2,1) = 1;
  
  // quadratic triangle
  Quad_Tri_En.Shape(3,3);
  Quad_Tri_En(0,0) = 1;
  Quad_Tri_En(0,1) = 4;
  Quad_Tri_En(0,2) = 2;
  Quad_Tri_En(1,0) = 2;
  Quad_Tri_En(1,1) = 5;
  Quad_Tri_En(1,2) = 3;
  Quad_Tri_En(2,0) = 3;
  Quad_Tri_En(2,1) = 6;
  Quad_Tri_En(2,2) = 1;
  
  Vel_Gn.Size(2);
  temp1 = 0.0;
  temp2 = 0.0;

  // Natural boundary condition
  // Fix this for Ibel
  Lbed = 2; // All velocity elements that have Nat BC are on edge 2

  Ln1 = Quad_Tri_En(Lbed,0) - 1;
  Ln2 = Quad_Tri_En(Lbed,2) - 1;

  NodLn1 = Nod(Gbe,Ln1) - 1;
  x1 = Glxy(NodLn1,0);
  y1 = Glxy(NodLn1,1);
  
  NodLn2 = Nod(Gbe,Ln2) - 1;
  x2 = Glxy(NodLn2,0);
  y2 = Glxy(NodLn2,1);
      
  h = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

  Sf.Size(3);
  
  for (int Ni = 0; Ni <= Ngp-1; Ni++)
  {
    Xi = GAUSPT(Ni,Ngp - 1);
    Const = GAUSWT(Ni,Ngp - 1) * h / 2.0; // Jacobian is h/2

    Sf1D(Xi,
	2,	// order of approximation is 2.
	Sf); 	// output

    s = Sf(1) * h / 2 + Sf(2) * h;
    // global coordinates of the barycentric Gauss points
    x = (1 - s / h) * x1 + (s / h) * x2;
    y = (1 - s / h) * y1 + (s / h) * y2;

    for (int j = 0; j <= 2; j++) // loop over node on this edge, only 3 nodes for quadratic element
    {
      jj = Nod(Gbe,Quad_Tri_En(Lbed,j)-1) - 1; // global node number of this local node

      if((jj + 1) % (2 * NX + 1) == 1) // Entrance
      {
	
	VelNatEntrance(x,y,Vel_Gn);

      }
      else if((jj + 1) % (2* NX + 1) == 0) // Exit
      {
	VelNatExit(x,y,Vel_Gn);
      }
      else if((jj + 1) <= 2 * NX + 1) // Lower plate
      {
	VelNatLowerPlate(x,y,Vel_Gn);
      }
      else // Upper plate
      {
	VelNatUpperPlate(x,y,Vel_Gn);
      }
      
      temp1 += Sf(j) * Vel_Gn(0) * Const;
      temp2 += Sf(j) * Vel_Gn(1) * Const;	
    }
  } // for NI
  
  Vel_Gn(0) = temp1;
  Vel_Gn(1) = temp2;
}