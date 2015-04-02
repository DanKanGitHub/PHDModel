#include "StrNodeDepartureFeet.h"

void StrNodeDepartureFeet(int myid,
			  int Nem,
			  int Vel_Npe,
			  int Str_Npe,
			  int VEL_FLAG,
			  int Vel_Nnm,
			  int NX,
			  int STRESS_FLAG,
			  double TIMESTEP,
			  E_SDM Vel_Glxy,
			  E_ISDM Vel_Nod,
			  E_ISDM Ele_Neigh,
			  E_ISDV All_Proc_Nodes,
			  E_ISDV My_Proc_Eles,
			  double *Vel_Old,
			  E_SDM & StrNodeDepartFootx,
			  E_SDM & StrNodeDepartFooty,
			  E_ISDM & StrNodeDepartElement,
		  E_ISDM Str_Nod,
		  E_SDM Str_Glxy,
		  E_SDM Str, 
		  E_SDM Str_Old)
{
  // Velocity variables
  E_SDM Vel_Elxy, Vel_Gdsf, El_Vel, Grad_Vel, Det_Inv_Grad_Vel;
  int Vel_Inod;
  
  // Pressure variables
  E_SDM Pre_Elxy;
  E_SDV Pre_Sf;
  
  // Stress variables
  E_SDM El_Dep_Str, El_Str, Str_Elxy;
  E_SDM Str_Gdsf;
  E_SDV It_Str, Str_Sf, Depart_Str; // iterative stress
  
  // Bio variables
  E_SDV El_Bio;
  double GaussPt_Bio;
  
  // Other variables
  E_SDM Gdsf, TempMat1, TempMat2, TempMat3, TempMat4, Dep_Gdsf;
  E_SDV Sf, Y_New_Xi_Eta, Ini_Foot, Depart_Foot, Dep_Sf;
  int Inod, New_Ele, Cur_Ele, JJ;//, Gauss_Pt_Num, jj, ii;
  double Xi, Eta, Coeff, Det_Grad_Vel, Trace_Grad_Vel;//Const1, Const2, ;
  double alpha1, alpha2, alpha3, x, y, x1, x2, x3, y1, y2, y3, Two_Area, DetJ;
  double RetardDivRelax, U_Zero, Effective_Pol_Vis, Effective_Den, ALPHA, BETA, P_Zero, Wi;
  double GaussPt_Bio_Weight, Effective_Vis;
  
  // Initialize
  // Velocity terms
  Vel_Elxy.Shape(Vel_Npe, 2);
  Vel_Gdsf.Shape(2, Vel_Npe);
  Det_Inv_Grad_Vel.Shape(2, 2);
  El_Vel.Shape(2, Vel_Npe);
  Sf.Size(Vel_Npe);  	// value of shape functions at (xi,eta)
  Grad_Vel.Shape(2, 2);

  // Stress terms
  El_Dep_Str.Shape(2,2);	// The Stress at the departure feet
  Str_Elxy.Shape(Str_Npe,2);
  It_Str.Size(Str_Npe);
  Str_Sf.Size(Str_Npe);
  Str_Gdsf.Shape(2,Str_Npe);
  Depart_Str.Size(Str_Npe);
  El_Str.Shape(4,Str_Npe);
  El_Bio.Size(Vel_Npe);
  
  // Other terms
  Y_New_Xi_Eta.Size(2);
  Ini_Foot.Size(2);
  Depart_Foot.Size(2);
  Dep_Sf.Size(Str_Npe);
  Dep_Gdsf.Shape(2,Str_Npe);

  for (int Ne = 0; Ne <= Nem-1; Ne++) 	// loop over all the elements
  {

    for (int i = 0; i <= Vel_Npe-1; i++) 	// get global coordinates of local nodes of element NE
    {
      Vel_Inod = Vel_Nod(Ne,i) - 1;		// Global node number (minus one for indexing) of local node.
      Vel_Elxy(i,0) = Vel_Glxy(Vel_Inod,0);	// x-coordinate
      Vel_Elxy(i,1) = Vel_Glxy(Vel_Inod,1);	// y-coordinate
    }

    // The first element
    if(Ne == 0)
    {
      for (int j = 0; j <= Str_Npe-1; j++)  // loop over all of the element nodes
      {
	Vel_Inod = Vel_Nod(Ne,j) - 1;

	x = Vel_Glxy(Vel_Inod, 0);
	y = Vel_Glxy(Vel_Inod, 1);

	Ini_Foot(0) = x;
	Ini_Foot(1) = y;

	DepartureFoot(VEL_FLAG, 
			Ne, 
			Vel_Nnm,
			TIMESTEP,
			Vel_Old, 
			Vel_Glxy, 
			Vel_Npe, 
			Vel_Nod, 
			Ele_Neigh,
			Ini_Foot,
			Depart_Foot,	// output
			New_Ele);	// output

	Cur_Ele = New_Ele;
	
	StrNodeDepartElement(Ne, j) = New_Ele;

	// Compute the stress at the departure foot
	for (int i = 0; i <= Str_Npe - 1; i++) 	// get global coordinates of local nodes of element NE
	{
	  Inod = Str_Nod(Cur_Ele-1,i) - 1;	// Global node number (minus one for C++ indexing) of local node.
	  Str_Elxy(i,0) = Str_Glxy(Inod,0);	// x-coordinate
	  Str_Elxy(i,1) = Str_Glxy(Inod,1);	// y-coordinate

	  El_Str(0,i) = Str_Old(Inod,0);
	  El_Str(1,i) = Str_Old(Inod,1);
	  El_Str(2,i) = El_Str(1,i); // The stress tensor is symmetric
	  El_Str(3,i) = Str_Old(Inod,2);
	}

	alpha1 = Str_Elxy(1,0) * Str_Elxy(2,1) - Str_Elxy(2,0) * Str_Elxy(1,1); // x2 * y3 - x3 * y2
	alpha2 = Str_Elxy(2,0) * Str_Elxy(0,1) - Str_Elxy(0,0) * Str_Elxy(2,1); // x3 * y1 - x1 * y3
	alpha3 = Str_Elxy(0,0) * Str_Elxy(1,1) - Str_Elxy(1,0) * Str_Elxy(0,1); // x1 * y2 - x2 * y1
	
	Two_Area = alpha1 + alpha2 + alpha3;
	
	Xi = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(2,0)) * (Str_Elxy(1,1) - Str_Elxy(2,1)) - 
			  (Depart_Foot(1) - Str_Elxy(2,1)) * (Str_Elxy(1,0) - Str_Elxy(2,0)));
	
	Eta = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(0,0)) * (Str_Elxy(2,1) - Str_Elxy(0,1)) + 
			  (Depart_Foot(1) - Str_Elxy(0,1)) * (Str_Elxy(0,0) - Str_Elxy(2,0)));
	
	StrNodeDepartFootx(Ne, j) = Xi;
	StrNodeDepartFooty(Ne, j) = Eta;

      } // for j
    }
    // minus 1 for C++ indexing Ele 2 given by Ne = 1.  Is Ne an odd numbered element on the first row.
    else if((Ne % 2 == 0) && (Ne <= 2 * NX - 2))
    {
      for (int j = 1; j <= Str_Npe - 1; j++)  // loop over nodes 2 through Str_Npe for the first row
      {
	Vel_Inod = Vel_Nod(Ne,j) - 1;

	x = Vel_Glxy(Vel_Inod, 0);
	y = Vel_Glxy(Vel_Inod, 1);

	Ini_Foot(0) = x;
	Ini_Foot(1) = y;

	DepartureFoot(VEL_FLAG, 
			Ne, 
			Vel_Nnm,
			TIMESTEP,
			Vel_Old, 
			Vel_Glxy, 
			Vel_Npe, 
			Vel_Nod, 
			Ele_Neigh,
			Ini_Foot,
			Depart_Foot,	// output
			New_Ele);		// output

	Cur_Ele = New_Ele;
	
	StrNodeDepartElement(Ne, j) = New_Ele;

	// Compute the stress at the departure foot
	for (int i = 0; i <= Str_Npe - 1; i++) 	// get global coordinates of local nodes of element NE
	{
	  Inod = Str_Nod(Cur_Ele-1,i) - 1;	// Global node number (minus one for C++ indexing) of local node.
	  Str_Elxy(i,0) = Str_Glxy(Inod,0);	// x-coordinate
	  Str_Elxy(i,1) = Str_Glxy(Inod,1);	// y-coordinate

	  El_Str(0,i) = Str_Old(Inod,0);
	  El_Str(1,i) = Str_Old(Inod,1);
	  El_Str(2,i) = El_Str(1,i); // The stress tensor is symmetric
	  El_Str(3,i) = Str_Old(Inod,2);
	}
	
	alpha1 = Str_Elxy(1,0) * Str_Elxy(2,1) - Str_Elxy(2,0) * Str_Elxy(1,1); // x2 * y3 - x3 * y2
	alpha2 = Str_Elxy(2,0) * Str_Elxy(0,1) - Str_Elxy(0,0) * Str_Elxy(2,1); // x3 * y1 - x1 * y3
	alpha3 = Str_Elxy(0,0) * Str_Elxy(1,1) - Str_Elxy(1,0) * Str_Elxy(0,1); // x1 * y2 - x2 * y1
	
	Two_Area = alpha1 + alpha2 + alpha3;
	
	Xi = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(2,0)) * (Str_Elxy(1,1) - Str_Elxy(2,1)) - 
			  (Depart_Foot(1) - Str_Elxy(2,1)) * (Str_Elxy(1,0) - Str_Elxy(2,0)));
	
	Eta = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(0,0)) * (Str_Elxy(2,1) - Str_Elxy(0,1)) + 
			  (Depart_Foot(1) - Str_Elxy(0,1)) * (Str_Elxy(0,0) - Str_Elxy(2,0)));
	
	StrNodeDepartFootx(Ne, j) = Xi;
	StrNodeDepartFooty(Ne, j) = Eta;

      } // for j
    }
      
    // If Ne is the first even element in a row or an odd element greater than those in the first row
    // The additional of the 2 * NX is to shift the first even elements value larger than 2 * NX - 1
    // Since this is modular arithmetic this shift has no impact on larger values of Ne.
    // or Ne is an odd element larger than the those in the first row.
    else if(((Ne + 1) % (2 * NX) == 2) || (Ne % 2 == 0))
    {
      JJ = 2; // Only solve at the 3rd node of each element
    
      Vel_Inod = Vel_Nod(Ne,JJ) - 1;

      x = Vel_Glxy(Vel_Inod, 0);
      y = Vel_Glxy(Vel_Inod, 1);

      Ini_Foot(0) = x;
      Ini_Foot(1) = y;

      DepartureFoot(VEL_FLAG, 
		      Ne, 
		      Vel_Nnm,
		      TIMESTEP,
		      Vel_Old, 
		      Vel_Glxy, 
		      Vel_Npe, 
		      Vel_Nod, 
		      Ele_Neigh,
		      Ini_Foot,
		      Depart_Foot,	// output
		      New_Ele);		// output

      Cur_Ele = New_Ele;
      
      StrNodeDepartElement(Ne, JJ) = New_Ele;

      // Compute the stress at the departure foot
      for (int i = 0; i <= Str_Npe - 1; i++) 	// get global coordinates of local nodes of element NE
      {
	Inod = Str_Nod(Cur_Ele-1,i) - 1;	// Global node number (minus one for C++ indexing) of local node.
	Str_Elxy(i,0) = Str_Glxy(Inod,0);	// x-coordinate
	Str_Elxy(i,1) = Str_Glxy(Inod,1);	// y-coordinate

	El_Str(0,i) = Str_Old(Inod,0);
	El_Str(1,i) = Str_Old(Inod,1);
	El_Str(2,i) = El_Str(1,i); // The stress tensor is symmetric
	El_Str(3,i) = Str_Old(Inod,2);
      }
      
      alpha1 = Str_Elxy(1,0) * Str_Elxy(2,1) - Str_Elxy(2,0) * Str_Elxy(1,1); // x2 * y3 - x3 * y2
      alpha2 = Str_Elxy(2,0) * Str_Elxy(0,1) - Str_Elxy(0,0) * Str_Elxy(2,1); // x3 * y1 - x1 * y3
      alpha3 = Str_Elxy(0,0) * Str_Elxy(1,1) - Str_Elxy(1,0) * Str_Elxy(0,1); // x1 * y2 - x2 * y1
      
      Two_Area = alpha1 + alpha2 + alpha3;
      
      Xi = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(2,0)) * (Str_Elxy(1,1) - Str_Elxy(2,1)) - 
			(Depart_Foot(1) - Str_Elxy(2,1)) * (Str_Elxy(1,0) - Str_Elxy(2,0)));
      
      Eta = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(0,0)) * (Str_Elxy(2,1) - Str_Elxy(0,1)) + 
			(Depart_Foot(1) - Str_Elxy(0,1)) * (Str_Elxy(0,0) - Str_Elxy(2,0)));
      
      StrNodeDepartFootx(Ne, JJ) = Xi;
      StrNodeDepartFooty(Ne, JJ) = Eta;

      if(STRESS_FLAG == 2)
      {
	for(int j = 4; j <= 5; j++)
	{
	  Vel_Inod = Vel_Nod(Ne,j) - 1;

	  x = Vel_Glxy(Vel_Inod, 0);
	  y = Vel_Glxy(Vel_Inod, 1);

	  Ini_Foot(0) = x;
	  Ini_Foot(1) = y;

	  DepartureFoot(VEL_FLAG, 
			  Ne, 
			  Vel_Nnm,
			  TIMESTEP,
			  Vel_Old, 
			  Vel_Glxy, 
			  Vel_Npe, 
			  Vel_Nod, 
			  Ele_Neigh,
			  Ini_Foot,
			  Depart_Foot,	// output
			  New_Ele);		// output

	  Cur_Ele = New_Ele;
	  
	  StrNodeDepartElement(Ne, j) = New_Ele;

	  // Compute the stress at the departure foot
	  for (int i = 0; i <= Str_Npe - 1; i++) 	// get global coordinates of local nodes of element NE
	  {
	    Inod = Str_Nod(Cur_Ele-1,i) - 1;	// Global node number (minus one for C++ indexing) of local node.
	    Str_Elxy(i,0) = Str_Glxy(Inod,0);	// x-coordinate
	    Str_Elxy(i,1) = Str_Glxy(Inod,1);	// y-coordinate

	    El_Str(0,i) = Str_Old(Inod,0);
	    El_Str(1,i) = Str_Old(Inod,1);
	    El_Str(2,i) = El_Str(1,i); // The stress tensor is symmetric
	    El_Str(3,i) = Str_Old(Inod,2);
	  }

	  alpha1 = Str_Elxy(1,0) * Str_Elxy(2,1) - Str_Elxy(2,0) * Str_Elxy(1,1); // x2 * y3 - x3 * y2
	  alpha2 = Str_Elxy(2,0) * Str_Elxy(0,1) - Str_Elxy(0,0) * Str_Elxy(2,1); // x3 * y1 - x1 * y3
	  alpha3 = Str_Elxy(0,0) * Str_Elxy(1,1) - Str_Elxy(1,0) * Str_Elxy(0,1); // x1 * y2 - x2 * y1
	  
	  Two_Area = alpha1 + alpha2 + alpha3;
	  
	  Xi = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(2,0)) * (Str_Elxy(1,1) - Str_Elxy(2,1)) - 
			    (Depart_Foot(1) - Str_Elxy(2,1)) * (Str_Elxy(1,0) - Str_Elxy(2,0)));
	  
	  Eta = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(0,0)) * (Str_Elxy(2,1) - Str_Elxy(0,1)) + 
			    (Depart_Foot(1) - Str_Elxy(0,1)) * (Str_Elxy(0,0) - Str_Elxy(2,0)));
	  
	  StrNodeDepartFootx(Ne, j) = Xi;
	  StrNodeDepartFooty(Ne, j) = Eta;
	}
      }
    }
    // If NX is not the first element in a row but is even
    else if(((Ne + 1) % (2 * NX) != 2) || NX % 2 == 1)
    {
      
      JJ = 4; // Only solve at the 5th node of each element
    
      Vel_Inod = Vel_Nod(Ne,JJ) - 1;

      x = Vel_Glxy(Vel_Inod, 0);
      y = Vel_Glxy(Vel_Inod, 1);

      Ini_Foot(0) = x;
      Ini_Foot(1) = y;

      DepartureFoot(VEL_FLAG, 
		      Ne, 
		      Vel_Nnm,
		      TIMESTEP,
		      Vel_Old, 
		      Vel_Glxy, 
		      Vel_Npe, 
		      Vel_Nod, 
		      Ele_Neigh,
		      Ini_Foot,
		      Depart_Foot,	// output
		      New_Ele);		// output

      Cur_Ele = New_Ele;
      
      StrNodeDepartElement(Ne, JJ) = New_Ele;

      // Compute the stress at the departure foot
      for (int i = 0; i <= Str_Npe - 1; i++) 	// get global coordinates of local nodes of element NE
      {
	Inod = Str_Nod(Cur_Ele-1,i) - 1;	// Global node number (minus one for C++ indexing) of local node.
	Str_Elxy(i,0) = Str_Glxy(Inod,0);	// x-coordinate
	Str_Elxy(i,1) = Str_Glxy(Inod,1);	// y-coordinate

	El_Str(0,i) = Str_Old(Inod,0);
	El_Str(1,i) = Str_Old(Inod,1);
	El_Str(2,i) = El_Str(1,i); // The stress tensor is symmetric
	El_Str(3,i) = Str_Old(Inod,2);
      }
      
      alpha1 = Str_Elxy(1,0) * Str_Elxy(2,1) - Str_Elxy(2,0) * Str_Elxy(1,1); // x2 * y3 - x3 * y2
      alpha2 = Str_Elxy(2,0) * Str_Elxy(0,1) - Str_Elxy(0,0) * Str_Elxy(2,1); // x3 * y1 - x1 * y3
      alpha3 = Str_Elxy(0,0) * Str_Elxy(1,1) - Str_Elxy(1,0) * Str_Elxy(0,1); // x1 * y2 - x2 * y1
      
      Two_Area = alpha1 + alpha2 + alpha3;
      
      Xi = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(2,0)) * (Str_Elxy(1,1) - Str_Elxy(2,1)) - 
			(Depart_Foot(1) - Str_Elxy(2,1)) * (Str_Elxy(1,0) - Str_Elxy(2,0)));
      
      Eta = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(0,0)) * (Str_Elxy(2,1) - Str_Elxy(0,1)) + 
			(Depart_Foot(1) - Str_Elxy(0,1)) * (Str_Elxy(0,0) - Str_Elxy(2,0)));
      
      StrNodeDepartFootx(Ne, JJ) = Xi;
      StrNodeDepartFooty(Ne, JJ) = Eta;

    }// end if
  } // for NE
}
// {
//   double Xi, Eta, alpha1, alpha2, alpha3, Two_Area, a1, a2, b1, b2, c1, c2;
//   
//   int Vel_Inod, New_Ele, Inod, jj, JJ;
//   
//   E_SDM Vel_Elxy;
//   
//   E_SDV Ini_Foot, Depart_Foot;
//   
//   Ini_Foot.Size(2);
//   Depart_Foot.Size(2);
//   Vel_Elxy.Shape(Vel_Npe, 2);
//   
//   for (int Ne = 0; Ne <= Nem - 1; Ne++) 		// loop over all the elements
//   {
//     for (int i = 0; i <= Vel_Npe - 1; i++) 	// get global coordinates of local nodes of element Ne
//     {
//       Vel_Inod = Vel_Nod(Ne,i) - 1;		// Global node number (minus one for indeXing) of local node.
//       Vel_Elxy(i,0) = Vel_Glxy(Vel_Inod,0);	// x-coordinate of the velcity
//       Vel_Elxy(i,1) = Vel_Glxy(Vel_Inod,1);	// y-coordinate
//     }
//     
//     for (int Npe = 0; Npe <= Str_Npe-1; Npe++) // loop over quadrature points
//     {
//       // This will work for both linear and quadratic stress since the first three nodes in each element are the same for both lin and quad
//       jj = Vel_Nod(Ne, Npe) - 1;
// 
//       // These are the global coordinates of the gauss point and are constant throughout this loop
//       Ini_Foot(0) = Vel_Glxy(jj,0); //Xi;
//       Ini_Foot(1) = Vel_Glxy(jj,1); //Eta;
// 
//       DepartureFoot(VEL_FLAG, 
// 		    Ne, 
// 		    Vel_Nnm,
// 		    TIMESTEP,
// 		    Vel_Old, 
// 		    Vel_Glxy, 
// 		    Vel_Npe, 
// 		    Vel_Nod, 
// 		    Ele_Neigh,
// 		    Ini_Foot,	// in x,y space
// 		    Depart_Foot,	// output in x,y space
// 		    New_Ele);	// output
// 
//       // Determine the element the Gauss point has moved into
//       if(New_Ele-1 != Ne)
//       {
// 	for (int k = 0; k <= Vel_Npe-1; k++) // get global coordinates of local nodes of element NE
// 	{
// 	  Inod = Vel_Nod(New_Ele-1,k) - 1;		// Global node number of local node.
// 	  Vel_Elxy(k,0) = Vel_Glxy(Inod,0);	// x-coordinate of te new element
// 	  Vel_Elxy(k,1) = Vel_Glxy(Inod,1);	// y-coordinate of the new element
// 	}
//       }
// 
//       // Convert the cartesian coordinates to barycentric coordinates
//       alpha1 = Vel_Elxy(1,0) * Vel_Elxy(2,1) - Vel_Elxy(2,0) * Vel_Elxy(1,1); // x2 * y3 - x3 * y2
//       alpha2 = Vel_Elxy(2,0) * Vel_Elxy(0,1) - Vel_Elxy(0,0) * Vel_Elxy(2,1); // x3 * y1 - x1 * y3
//       alpha3 = Vel_Elxy(0,0) * Vel_Elxy(1,1) - Vel_Elxy(1,0) * Vel_Elxy(0,1); // x1 * y2 - x2 * y1
//       
//       Two_Area = alpha1 + alpha2 + alpha3;
//       
//       StrNodeDepartFootx(Ne, Npe) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(2,0)) * (Vel_Elxy(1,1) - Vel_Elxy(2,1)) - 
// 			(Depart_Foot(1) - Vel_Elxy(2,1)) * (Vel_Elxy(1,0) - Vel_Elxy(2,0)));
//       
//       StrNodeDepartFooty(Ne, Npe) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(0,0)) * (Vel_Elxy(2,1) - Vel_Elxy(0,1)) + 
// 			(Depart_Foot(1) - Vel_Elxy(0,1)) * (Vel_Elxy(0,0) - Vel_Elxy(2,0)));
//       
//       StrNodeDepartElement(Ne, Npe) = New_Ele;
//     } // for Npe
    // The first element
//     if(Ne == 0)
//     {
//       for (int j = 0; j <= Str_Npe-1; j++)  // loop over all of the element nodes
//       {
// 	
// 	// This will work for both linear and quadratic stress since the first three nodes in each element are the same for both lin and quad
// 	jj = Vel_Nod(Ne, j) - 1;
// 
// 	// These are the global coordinates of the gauss point and are constant throughout this loop
// 	Ini_Foot(0) = Vel_Glxy(jj,0); //Xi;
// 	Ini_Foot(1) = Vel_Glxy(jj,1); //Eta;
// 
// 	DepartureFoot(VEL_FLAG, 
// 		      Ne, 
// 		      Vel_Nnm,
// 		      TIMESTEP,
// 		      Vel_Old, 
// 		      Vel_Glxy, 
// 		      Vel_Npe, 
// 		      Vel_Nod, 
// 		      Ele_Neigh,
// 		      Ini_Foot,	// in x,y space
// 		      Depart_Foot,	// output in x,y space
// 		      New_Ele);	// output
// 
// 	// Determine the element the Gauss point has moved into
// // 	if(New_Ele-1 != Ne)
// // 	{
// 	  for (int k = 0; k <= Vel_Npe-1; k++) // get global coordinates of local nodes of element NE
// 	  {
// 	    Inod = Vel_Nod(New_Ele-1,k) - 1;		// Global node number of local node.
// 	    Vel_Elxy(k,0) = Vel_Glxy(Inod,0);	// x-coordinate of te new element
// 	    Vel_Elxy(k,1) = Vel_Glxy(Inod,1);	// y-coordinate of the new element
// 	  }
// // 	}
// 
// 	// Convert the cartesian coordinates to barycentric coordinates
// 	alpha1 = Vel_Elxy(1,0) * Vel_Elxy(2,1) - Vel_Elxy(2,0) * Vel_Elxy(1,1); // x2 * y3 - x3 * y2
// 	alpha2 = Vel_Elxy(2,0) * Vel_Elxy(0,1) - Vel_Elxy(0,0) * Vel_Elxy(2,1); // x3 * y1 - x1 * y3
// 	alpha3 = Vel_Elxy(0,0) * Vel_Elxy(1,1) - Vel_Elxy(1,0) * Vel_Elxy(0,1); // x1 * y2 - x2 * y1
// 	
// 	Two_Area = alpha1 + alpha2 + alpha3;
// 	
// 	StrNodeDepartFootx(Ne, j) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(2,0)) * (Vel_Elxy(1,1) - Vel_Elxy(2,1)) - 
// 			  (Depart_Foot(1) - Vel_Elxy(2,1)) * (Vel_Elxy(1,0) - Vel_Elxy(2,0)));
// 	
// 	StrNodeDepartFooty(Ne, j) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(0,0)) * (Vel_Elxy(2,1) - Vel_Elxy(0,1)) + 
// 			  (Depart_Foot(1) - Vel_Elxy(0,1)) * (Vel_Elxy(0,0) - Vel_Elxy(2,0)));
// 	
// 	StrNodeDepartElement(Ne, j) = New_Ele;
// 
//       } // for j
//     } // if Ne
//     // minus 1 for C++ indexing Ele 2 given by Ne = 1.  Is Ne an odd numbered element on the first row.
//     else if((Ne % 2 == 0) && (Ne <= 2 * NX - 2))
//     {
//       for (int j = 1; j <= Str_Npe - 1; j++)  // loop over nodes 2 through Str_Npe for the first row
//       {
// 	
// 	// This will work for both linear and quadratic stress since the first three nodes in each element are the same for both lin and quad
// 	jj = Vel_Nod(Ne, j) - 1;
// 
// 	// These are the global coordinates of the gauss point and are constant throughout this loop
// 	Ini_Foot(0) = Vel_Glxy(jj,0); //Xi;
// 	Ini_Foot(1) = Vel_Glxy(jj,1); //Eta;
// 
// 	DepartureFoot(VEL_FLAG, 
// 		      Ne, 
// 		      Vel_Nnm,
// 		      TIMESTEP,
// 		      Vel_Old, 
// 		      Vel_Glxy, 
// 		      Vel_Npe, 
// 		      Vel_Nod, 
// 		      Ele_Neigh,
// 		      Ini_Foot,	// in x,y space
// 		      Depart_Foot,	// output in x,y space
// 		      New_Ele);	// output
// 
// 	// Determine the element the Gauss point has moved into
// // 	if(New_Ele-1 != Ne)
// // 	{
// 	  for (int k = 0; k <= Vel_Npe-1; k++) // get global coordinates of local nodes of element NE
// 	  {
// 	    Inod = Vel_Nod(New_Ele-1,k) - 1;		// Global node number of local node.
// 	    Vel_Elxy(k,0) = Vel_Glxy(Inod,0);	// x-coordinate of te new element
// 	    Vel_Elxy(k,1) = Vel_Glxy(Inod,1);	// y-coordinate of the new element
// 	  }
// // 	}
// 
// 	// Convert the cartesian coordinates to barycentric coordinates
// 	alpha1 = Vel_Elxy(1,0) * Vel_Elxy(2,1) - Vel_Elxy(2,0) * Vel_Elxy(1,1); // x2 * y3 - x3 * y2
// 	alpha2 = Vel_Elxy(2,0) * Vel_Elxy(0,1) - Vel_Elxy(0,0) * Vel_Elxy(2,1); // x3 * y1 - x1 * y3
// 	alpha3 = Vel_Elxy(0,0) * Vel_Elxy(1,1) - Vel_Elxy(1,0) * Vel_Elxy(0,1); // x1 * y2 - x2 * y1
// 	
// 	Two_Area = alpha1 + alpha2 + alpha3;
// 	
// 	StrNodeDepartFootx(Ne, j) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(2,0)) * (Vel_Elxy(1,1) - Vel_Elxy(2,1)) - 
// 			  (Depart_Foot(1) - Vel_Elxy(2,1)) * (Vel_Elxy(1,0) - Vel_Elxy(2,0)));
// 	
// 	StrNodeDepartFooty(Ne, j) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(0,0)) * (Vel_Elxy(2,1) - Vel_Elxy(0,1)) + 
// 			  (Depart_Foot(1) - Vel_Elxy(0,1)) * (Vel_Elxy(0,0) - Vel_Elxy(2,0)));
// 	
// 	StrNodeDepartElement(Ne, j) = New_Ele;
// 
//       } // for j
//     }
//       
//     // If Ne is the first even element in a row or an odd element greater than those in the first row
//     // The additional of the 2 * NX is to shift the first even elements value larger than 2 * NX - 1
//     // Since this is modular arithmetic this shift has no impact on larger values of Ne.
//     // or Ne is an odd element larger than the those in the first row.
//     else if(((Ne + 1) % (2 * NX) == 2) || (Ne % 2 == 0))
//     {
// 
//       JJ = 2; // Only solve at the 3rd node of each element
// 
//       // These are the global coordinates of the gauss point and are constant throughout this loop
//       Ini_Foot(0) = Vel_Glxy(JJ,0); //Xi;
//       Ini_Foot(1) = Vel_Glxy(JJ,1); //Eta;
// 
//       DepartureFoot(VEL_FLAG, 
// 		    Ne, 
// 		    Vel_Nnm,
// 		    TIMESTEP,
// 		    Vel_Old, 
// 		    Vel_Glxy, 
// 		    Vel_Npe, 
// 		    Vel_Nod, 
// 		    Ele_Neigh,
// 		    Ini_Foot,	// in x,y space
// 		    Depart_Foot,	// output in x,y space
// 		    New_Ele);	// output
// 
//       // Determine the element the Gauss point has moved into
// //       if(New_Ele-1 != Ne)
// //       {
// 	for (int k = 0; k <= Vel_Npe-1; k++) // get global coordinates of local nodes of element NE
// 	{
// 	  Inod = Vel_Nod(New_Ele-1,k) - 1;		// Global node number of local node.
// 	  Vel_Elxy(k,0) = Vel_Glxy(Inod,0);	// x-coordinate of te new element
// 	  Vel_Elxy(k,1) = Vel_Glxy(Inod,1);	// y-coordinate of the new element
// 	}
// //       }
// 
//       // Convert the cartesian coordinates to barycentric coordinates
//       alpha1 = Vel_Elxy(1,0) * Vel_Elxy(2,1) - Vel_Elxy(2,0) * Vel_Elxy(1,1); // x2 * y3 - x3 * y2
//       alpha2 = Vel_Elxy(2,0) * Vel_Elxy(0,1) - Vel_Elxy(0,0) * Vel_Elxy(2,1); // x3 * y1 - x1 * y3
//       alpha3 = Vel_Elxy(0,0) * Vel_Elxy(1,1) - Vel_Elxy(1,0) * Vel_Elxy(0,1); // x1 * y2 - x2 * y1
//       
//       Two_Area = alpha1 + alpha2 + alpha3;
//       
//       StrNodeDepartFootx(Ne, JJ) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(2,0)) * (Vel_Elxy(1,1) - Vel_Elxy(2,1)) - 
// 			(Depart_Foot(1) - Vel_Elxy(2,1)) * (Vel_Elxy(1,0) - Vel_Elxy(2,0)));
//       
//       StrNodeDepartFooty(Ne, JJ) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(0,0)) * (Vel_Elxy(2,1) - Vel_Elxy(0,1)) + 
// 			(Depart_Foot(1) - Vel_Elxy(0,1)) * (Vel_Elxy(0,0) - Vel_Elxy(2,0)));
//       
//       StrNodeDepartElement(Ne, JJ) = New_Ele;
// 
//       if(STRESS_FLAG == 2)
//       {
// 	for(int j = 4; j <= 5; j++)
// 	{
// 	  
// 	  // This will work for both linear and quadratic stress since the first three nodes in each element are the same for both lin and quad
// 	  jj = Vel_Nod(Ne, j) - 1;
// 
// 	  // These are the global coordinates of the gauss point and are constant throughout this loop
// 	  Ini_Foot(0) = Vel_Glxy(jj,0); //Xi;
// 	  Ini_Foot(1) = Vel_Glxy(jj,1); //Eta;
// 
// 	  DepartureFoot(VEL_FLAG, 
// 			Ne, 
// 			Vel_Nnm,
// 			TIMESTEP,
// 			Vel_Old, 
// 			Vel_Glxy, 
// 			Vel_Npe, 
// 			Vel_Nod, 
// 			Ele_Neigh,
// 			Ini_Foot,	// in x,y space
// 			Depart_Foot,	// output in x,y space
// 			New_Ele);	// output
// 
// 	  // Determine the element the Gauss point has moved into
// // 	  if(New_Ele-1 != Ne)
// // 	  {
// 	    for (int k = 0; k <= Vel_Npe-1; k++) // get global coordinates of local nodes of element NE
// 	    {
// 	      Inod = Vel_Nod(New_Ele-1,k) - 1;		// Global node number of local node.
// 	      Vel_Elxy(k,0) = Vel_Glxy(Inod,0);	// x-coordinate of te new element
// 	      Vel_Elxy(k,1) = Vel_Glxy(Inod,1);	// y-coordinate of the new element
// 	    }
// // 	  }
// 
// 	  // Convert the cartesian coordinates to barycentric coordinates
// 	  alpha1 = Vel_Elxy(1,0) * Vel_Elxy(2,1) - Vel_Elxy(2,0) * Vel_Elxy(1,1); // x2 * y3 - x3 * y2
// 	  alpha2 = Vel_Elxy(2,0) * Vel_Elxy(0,1) - Vel_Elxy(0,0) * Vel_Elxy(2,1); // x3 * y1 - x1 * y3
// 	  alpha3 = Vel_Elxy(0,0) * Vel_Elxy(1,1) - Vel_Elxy(1,0) * Vel_Elxy(0,1); // x1 * y2 - x2 * y1
// 	  
// 	  Two_Area = alpha1 + alpha2 + alpha3;
// 	  
// 	  StrNodeDepartFootx(Ne, j) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(2,0)) * (Vel_Elxy(1,1) - Vel_Elxy(2,1)) - 
// 			    (Depart_Foot(1) - Vel_Elxy(2,1)) * (Vel_Elxy(1,0) - Vel_Elxy(2,0)));
// 	  
// 	  StrNodeDepartFooty(Ne, j) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(0,0)) * (Vel_Elxy(2,1) - Vel_Elxy(0,1)) + 
// 			    (Depart_Foot(1) - Vel_Elxy(0,1)) * (Vel_Elxy(0,0) - Vel_Elxy(2,0)));
// 	  
// 	  StrNodeDepartElement(Ne, j) = New_Ele;
// 
// 	} // for
//       } // if Stressflag
//     }
//     // If NX is not the first element in a row but is even
//     else if(((Ne + 1) % (2 * NX) != 2) || NX % 2 == 1)
//     {
// 
//       JJ = 4; // Only solve at the 5th node of each element
// 
//       // These are the global coordinates of the gauss point and are constant throughout this loop
//       Ini_Foot(0) = Vel_Glxy(JJ,0); //Xi;
//       Ini_Foot(1) = Vel_Glxy(JJ,1); //Eta;
// 
//       DepartureFoot(VEL_FLAG, 
// 		    Ne, 
// 		    Vel_Nnm,
// 		    TIMESTEP,
// 		    Vel_Old, 
// 		    Vel_Glxy, 
// 		    Vel_Npe, 
// 		    Vel_Nod, 
// 		    Ele_Neigh,
// 		    Ini_Foot,	// in x,y space
// 		    Depart_Foot,	// output in x,y space
// 		    New_Ele);	// output
// 
//       // Determine the element the Gauss point has moved into
// //       if(New_Ele-1 != Ne)
// //       {
// 	for (int k = 0; k <= Vel_Npe-1; k++) // get global coordinates of local nodes of element NE
// 	{
// 	  Inod = Vel_Nod(New_Ele-1,k) - 1;		// Global node number of local node.
// 	  Vel_Elxy(k,0) = Vel_Glxy(Inod,0);	// x-coordinate of te new element
// 	  Vel_Elxy(k,1) = Vel_Glxy(Inod,1);	// y-coordinate of the new element
// 	}
// //       }
// 
//       // Convert the cartesian coordinates to barycentric coordinates
//       alpha1 = Vel_Elxy(1,0) * Vel_Elxy(2,1) - Vel_Elxy(2,0) * Vel_Elxy(1,1); // x2 * y3 - x3 * y2
//       alpha2 = Vel_Elxy(2,0) * Vel_Elxy(0,1) - Vel_Elxy(0,0) * Vel_Elxy(2,1); // x3 * y1 - x1 * y3
//       alpha3 = Vel_Elxy(0,0) * Vel_Elxy(1,1) - Vel_Elxy(1,0) * Vel_Elxy(0,1); // x1 * y2 - x2 * y1
//       
//       Two_Area = alpha1 + alpha2 + alpha3;
//       
//       StrNodeDepartFootx(Ne, JJ) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(2,0)) * (Vel_Elxy(1,1) - Vel_Elxy(2,1)) - 
// 			(Depart_Foot(1) - Vel_Elxy(2,1)) * (Vel_Elxy(1,0) - Vel_Elxy(2,0)));
//       
//       StrNodeDepartFooty(Ne, JJ) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(0,0)) * (Vel_Elxy(2,1) - Vel_Elxy(0,1)) + 
// 			(Depart_Foot(1) - Vel_Elxy(0,1)) * (Vel_Elxy(0,0) - Vel_Elxy(2,0)));
//       
//       StrNodeDepartElement(Ne, JJ) = New_Ele;
// 
//     } // if Ne
//   } // for Ne
// }
