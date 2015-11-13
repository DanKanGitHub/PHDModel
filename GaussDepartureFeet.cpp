// The output of this subroutine is in barysentric coordinates
#include "GaussDepartureFeet.h"

void GaussDepartureFeet(int myid,
			int Nem,
			int N_TRI_QUAD,
			int Vel_Npe,
			int VEL_FLAG,
			int Vel_Nnm,
			int NX,
			double TIMESTEP,
			double dx,
			double dy,
			E_SDM Vel_Glxy,
			E_SDM Tri_Quad_Pt,
			E_ISDM Vel_Nod,
			E_ISDM Ele_Neigh,
			E_ISDV My_Proc_Eles,
			double *Vel_Old,
			E_SDM & GaussDepartFootx,
			E_SDM & GaussDepartFooty,
			E_ISDM & GaussDepartElement)
{
  double Xi, Eta, alpha1, alpha2, alpha3, Two_Area, a1, a2, b1, b2, c1, c2;
  
  int Vel_Inod, New_Ele, Inod;
  
  E_SDM Vel_Elxy;
  
  E_SDV Ini_Foot, Depart_Foot;
  
  Ini_Foot.Size(2);
  Depart_Foot.Size(2);
  Vel_Elxy.Shape(Vel_Npe, 2);
  
  for (int Ne = 0; Ne <= Nem - 1; Ne++) 		// loop over all the elements
  {
    if(My_Proc_Eles(Ne) == myid)
    {
      for (int i = 0; i <= Vel_Npe - 1; i++) 	// get global coordinates of local nodes of element Ne
      {
	Vel_Inod = Vel_Nod(Ne,i) - 1;		// Global node number (minus one for indeXing) of local node.
	Vel_Elxy(i,0) = Vel_Glxy(Vel_Inod,0);	// x-coordinate of the velcity
	Vel_Elxy(i,1) = Vel_Glxy(Vel_Inod,1);	// y-coordinate
      }
      
      for (int Ni = 0; Ni <= N_TRI_QUAD-1; Ni++) // loop over quadrature points
      {
	Xi  = Tri_Quad_Pt(Ni,0);
	Eta = Tri_Quad_Pt(Ni,1);

	// Convert the gauss point into cartesian coordinates for the associated element.
	// In general: 2 * Area = x2 * y3 - x3 * y2 + x3 * y1 - x1 * y3 + x1 * y2 - x2 * y1
	alpha1 = Vel_Elxy(1,0) * Vel_Elxy(2,1) - Vel_Elxy(2,0) * Vel_Elxy(1,1); // x2 * y3 - x3 * y2
	alpha2 = Vel_Elxy(2,0) * Vel_Elxy(0,1) - Vel_Elxy(0,0) * Vel_Elxy(2,1); // x3 * y1 - x1 * y3
	alpha3 = Vel_Elxy(0,0) * Vel_Elxy(1,1) - Vel_Elxy(1,0) * Vel_Elxy(0,1); // x1 * y2 - x2 * y1
	
	Two_Area = alpha1 + alpha2 + alpha3;

	a1 = Two_Area * Xi + Vel_Elxy(2,0) * (Vel_Elxy(1,1) - Vel_Elxy(2,1)) - Vel_Elxy(2,1) * (Vel_Elxy(1,0) - Vel_Elxy(2,0));

	b1 = Vel_Elxy(1,1) - Vel_Elxy(2,1);

	c1 = Vel_Elxy(2,0) - Vel_Elxy(1,0);

	a2 = Two_Area * Eta + Vel_Elxy(0,0) * (Vel_Elxy(2,1) - Vel_Elxy(0,1)) + Vel_Elxy(0,1) * (Vel_Elxy(0,0) - Vel_Elxy(2,0));

	b2 = Vel_Elxy(2,1) - Vel_Elxy(0,1);

	c2 = Vel_Elxy(0,0) - Vel_Elxy(2,0);

	// These are the global coordinates of the gauss point and are constant throughout this loop
	Ini_Foot(0) = 1.0 / (b1 * c2 - c1 * b2) * (c2 * a1 - c1 * a2);

	Ini_Foot(1) = 1.0 / (b1 * c2 - c1 * b2) * (b1 * a2 - b2 * a1);
	
// 	std::cout << "Ini_Foot(0) = " << Ini_Foot(0) << endl;
// 	std::cout << "Ini_Foot(1) = " << Ini_Foot(1) << endl;
// 	std::cout << "Ne = " << Ne << endl;

	DepartureFoot(VEL_FLAG, 
		      Ne, 
		      Vel_Nnm,
		      NX,
		      TIMESTEP,
		      dx,
		      dy,
		      Vel_Old, 
		      Vel_Glxy, 
		      Vel_Npe, 
		      Vel_Nod, 
		      Ele_Neigh,
		      Ini_Foot,		// in x,y space
		      Depart_Foot,	// output in x,y space
		      New_Ele);		// output
	
// 	int QWERTY;
// 	std::cin >> QWERTY;

	// Determine the element the Gauss point has moved into
	if(New_Ele-1 != Ne)
	{
	  for (int k = 0; k <= Vel_Npe-1; k++) // get global coordinates of local nodes of element NE
	  {
	    Inod = Vel_Nod(New_Ele-1,k) - 1;		// Global node number of local node.
	    Vel_Elxy(k,0) = Vel_Glxy(Inod,0);	// x-coordinate of te new element
	    Vel_Elxy(k,1) = Vel_Glxy(Inod,1);	// y-coordinate of the new element
	  }
	}

	// Convert the cartesian coordinates to barycentric coordinates
	alpha1 = Vel_Elxy(1,0) * Vel_Elxy(2,1) - Vel_Elxy(2,0) * Vel_Elxy(1,1); // x2 * y3 - x3 * y2
	alpha2 = Vel_Elxy(2,0) * Vel_Elxy(0,1) - Vel_Elxy(0,0) * Vel_Elxy(2,1); // x3 * y1 - x1 * y3
	alpha3 = Vel_Elxy(0,0) * Vel_Elxy(1,1) - Vel_Elxy(1,0) * Vel_Elxy(0,1); // x1 * y2 - x2 * y1
	
	Two_Area = alpha1 + alpha2 + alpha3;
	
	GaussDepartFootx(Ne, Ni) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(2,0)) * (Vel_Elxy(1,1) - Vel_Elxy(2,1)) - 
			 (Depart_Foot(1) - Vel_Elxy(2,1)) * (Vel_Elxy(1,0) - Vel_Elxy(2,0)));
	
	GaussDepartFooty(Ne, Ni) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(0,0)) * (Vel_Elxy(2,1) - Vel_Elxy(0,1)) + 
			  (Depart_Foot(1) - Vel_Elxy(0,1)) * (Vel_Elxy(0,0) - Vel_Elxy(2,0)));
	
	GaussDepartElement(Ne, Ni) = New_Ele;

      } // for Ni
    } // If My_Proc_eles
  } // for Ne
}