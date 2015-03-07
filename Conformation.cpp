// conformation tensor solver.  Uses linear elements only
// Can reduce the number of Gauss points for this program and have the same accuracy

//   Wi = RELAX_TIME / T_ZERO;	// Fix This !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//   Alpha = 1.0 / Wi;	
// Effective_Pol_Vis = (1 - RETARD_TIME / RELAX_TIME) * POL_NEW_VIS / (T_ZERO * P_Zero); // Non-Dim correction is / (T_ZERO * P_Zero) // not used in any subroutines
// Beta = Effective_Pol_Vis / (Wi * Wi);	

#include <cmath>
#include "Shape2d.h"
#include "DepartureFoot.h"
#include "Conformation.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_IntSerialDenseMatrix E_ISDM;
typedef Epetra_SerialDenseVector E_SDV;

void Conformation(double TIMESTEP, 
		  double Sol_Vis,
		  double Poly_Vis,
		  int NX,
		  double RETARD_TIME,
		  double Sol_Density,
		  double Poly_Density,
		  double T_ZERO,
		  double L_ZERO,
		  int Vel_Npe, 
		  int Pre_Npe, 
		  int Nem, 
		  int Vel_Nnm,
		  int N_TRI_QUAD,
		  int Pre_Nnm,
		  E_ISDM Vel_Nod, 
		  E_ISDM Pre_Nod, 
		  E_SDM Vel_Glxy, 
		  E_SDM Pre_Glxy, 
		  E_SDM Tri_Quad_Pt,
		  E_SDV Tri_Quad_Wt, 
		  E_ISDM Ele_Neigh,
		  int VEL_FLAG, 
		  int PRE_FLAG, 
		  double *Vel, 
		  double *Vel_Old,
		  double *Bio,
		  E_SDM Str, 
		  E_SDM Str_Old,
		  E_SDM & Str_New)
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
  
  // Pressure terms
  Pre_Elxy.Shape(Pre_Npe, 2);
  Pre_Sf.Size(Pre_Npe);
  Gdsf.Shape(2, Pre_Npe); 	// derivatives w.r.t. global cooridinates
  
  // Stress terms
  El_Dep_Str.Shape(2,2);	// The Stress at the departure feet
  Str_Elxy.Shape(Pre_Npe,2);
  It_Str.Size(3);
  Str_Sf.Size(Pre_Npe);
  Str_Gdsf.Shape(2,Pre_Npe);
  Depart_Str.Size(3);
  El_Str.Shape(4,Pre_Npe);
  El_Bio.Size(Vel_Npe);
  
  // Other terms
  TempMat1.Shape(2,2);
  TempMat2.Shape(2,2);
  TempMat3.Shape(2,2);
  TempMat4.Shape(2,2);
  Y_New_Xi_Eta.Size(2);
  Ini_Foot.Size(2);
  Depart_Foot.Size(2);
  Dep_Sf.Size(Pre_Npe);
  Dep_Gdsf.Shape(2,Pre_Npe);

  U_Zero = L_ZERO / T_ZERO;
  
  for (int Ne = 0; Ne <= Nem-1; Ne++) 	// loop over all the elements
  {
    for (int i = 0; i <= Pre_Npe-1; i++) // get global coordinates of local nodes of element NE
    {
      Inod = Pre_Nod(Ne,i) - 1;		// Global node number (minus one for indexing) of local node.
      Pre_Elxy(i,0) = Pre_Glxy(Inod,0);  // x-coordinate
      Pre_Elxy(i,1) = Pre_Glxy(Inod,1);  // y-coordinate
    }
      
    for (int i = 0; i <= Vel_Npe-1; i++) 	// get global coordinates of local nodes of element NE
    {
      Vel_Inod = Vel_Nod(Ne,i) - 1;		// Global node number (minus one for indexing) of local node.
      Vel_Elxy(i,0) = Vel_Glxy(Vel_Inod,0);	// x-coordinate
      Vel_Elxy(i,1) = Vel_Glxy(Vel_Inod,1);	// y-coordinate

      El_Vel(0,i) = Vel[Vel_Inod];
      El_Vel(1,i) = Vel[Vel_Inod + Vel_Nnm];
      
      El_Bio(i) = Bio[Vel_Inod];
    }

    if(Ne == 0)
    {
      for (int j = 0; j <= Pre_Npe-1; j++)  // loop over element nodes
      {
	Vel_Inod = Vel_Nod(Ne,j) - 1;

	alpha1 = Vel_Elxy(1,0) * Vel_Elxy(2,1) - Vel_Elxy(2,0) * Vel_Elxy(1,1); // x2 * y3 - x3 * y2
	alpha2 = Vel_Elxy(2,0) * Vel_Elxy(0,1) - Vel_Elxy(0,0) * Vel_Elxy(2,1); // x3 * y1 - x1 * y3
	alpha3 = Vel_Elxy(0,0) * Vel_Elxy(1,1) - Vel_Elxy(1,0) * Vel_Elxy(0,1); // x1 * y2 - x2 * y1
	
	// In general: 2 * Area = x2 * y3 - x3 * y2 + x3 * y1 - x1 * y3 + x1 * y2 - x2 * y1
	Two_Area = alpha1 + alpha2 + alpha3;
	x = Vel_Glxy(Vel_Inod, 0);
	y = Vel_Glxy(Vel_Inod, 1);
	x1 = Pre_Elxy(0,0);
	x2 = Pre_Elxy(1,0);
	x3 = Pre_Elxy(2,0);
	y1 = Pre_Elxy(0,1);
	y2 = Pre_Elxy(1,1);
	y3 = Pre_Elxy(2,1);
	
	Xi  = 1.0 / Two_Area * ((x - x3) * (y2 - y3) - (y - y3) * (x2 - x3));
	Eta = 1.0 / Two_Area * ((x - x1) * (y3 - y1) + (y - y1) * (x1 - x3));

	Shape2d(Xi, 
		Eta, 
		Vel_Elxy, 
		Vel_Npe, 
		VEL_FLAG, 
		Sf, 
		Vel_Gdsf, 
		DetJ);

	Shape2d(Xi, 
		Eta, 
		Pre_Elxy, 
		Pre_Npe, 
		PRE_FLAG, 
		Pre_Sf, 
		Gdsf, 
		DetJ);

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

	Cur_Ele = Ne + 1;
	
	// Compute the stress at the departure foot
	for (int i = 0; i <= Pre_Npe - 1; i++) 	// get global coordinates of local nodes of element NE
	{
	  Inod = Pre_Nod(Cur_Ele-1,i) - 1;	// Global node number (minus one for C++ indexing) of local node.
	  Str_Elxy(i,0) = Pre_Glxy(Inod,0);	// x-coordinate
	  Str_Elxy(i,1) = Pre_Glxy(Inod,1);	// y-coordinate

	  El_Str(0,i) = Str_Old(Inod,0);
	  El_Str(1,i) = Str_Old(Inod,1);
	  El_Str(2,i) = El_Str(1,i); // The stress tensor is symmetric
	  El_Str(3,i) = Str_Old(Inod,2);

	}

	alpha1 = Str_Elxy(1,0) * Str_Elxy(2,1) - Str_Elxy(2,0) * Str_Elxy(1,1); // x2 * y3 - x3 * y2
	alpha2 = Str_Elxy(2,0) * Str_Elxy(0,1) - Str_Elxy(0,0) * Str_Elxy(2,1); // x3 * y1 - x1 * y3
	alpha3 = Str_Elxy(0,0) * Str_Elxy(1,1) - Str_Elxy(1,0) * Str_Elxy(0,1); // x1 * y2 - x2 * y1
	
	Two_Area = alpha1 + alpha2 + alpha3;
	
	Y_New_Xi_Eta(0) = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(2,0)) * (Str_Elxy(1,1) - Str_Elxy(2,1)) - 
			  (Depart_Foot(1) - Str_Elxy(2,1)) * (Str_Elxy(1,0) - Str_Elxy(2,0)));
	
	Y_New_Xi_Eta(1) = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(0,0)) * (Str_Elxy(2,1) - Str_Elxy(0,1)) + 
			  (Depart_Foot(1) - Str_Elxy(0,1)) * (Str_Elxy(0,0) - Str_Elxy(2,0)));
	
	//The point may have not left the element but it did move so recompute everything at the new position
	Xi = Y_New_Xi_Eta(0);
	Eta = Y_New_Xi_Eta(1);

	Shape2d(Xi, 
		Eta, 
		Str_Elxy, 
		Pre_Npe, 	// Both stress and pressure are linear
		PRE_FLAG, 
		Dep_Sf, 	// output
		Dep_Gdsf, // output
		DetJ);	// output

	for(int i = 0; i <= 2; i++)
	{
	  It_Str(i) = 0.0;
	}

	for(int k = 0; k <= Pre_Npe - 1; k++)
	{
	  It_Str(0) += El_Str(0,k) * Dep_Sf(k);
	  It_Str(1) += El_Str(1,k) * Dep_Sf(k);
	  It_Str(2) += El_Str(3,k) * Dep_Sf(k); // Since the Str is symmetric
	}
	
	Depart_Str(0) = It_Str(0);
	Depart_Str(1) = It_Str(1);
	Depart_Str(2) = It_Str(2);
	
	// This is not the element departure stress but the Stress at the particular departed Guass point
	El_Dep_Str(0,0) = Depart_Str(0);
	El_Dep_Str(0,1) = Depart_Str(1);
	El_Dep_Str(1,0) = El_Dep_Str(0,1); // The stress tensor is symmetric
	El_Dep_Str(1,1) = Depart_Str(2);

	GaussPt_Bio = 0.0;
	
	for(int k = 0; k <= Vel_Npe - 1; k++)
	{
	  GaussPt_Bio += El_Bio(k) * Sf(k);
	}
	
	GaussPt_Bio_Weight = BioWeightFunc(GaussPt_Bio);
	
	RetardDivRelax = RetardationDividedByRelaxation(GaussPt_Bio_Weight);

	Effective_Den = (1 - GaussPt_Bio_Weight) * Sol_Density + GaussPt_Bio_Weight * Poly_Density;
	
	P_Zero = Effective_Den * U_Zero * U_Zero;
	
	Effective_Vis = (Sol_Vis * (1 - GaussPt_Bio_Weight) + Poly_Vis * GaussPt_Bio_Weight) / (T_ZERO * P_Zero);
	
	Effective_Pol_Vis = (1 - RetardDivRelax) * Effective_Vis;
	
	Wi = RETARD_TIME / RetardDivRelax; // RETARD_TIME / RetardDivRelax = RELAX_TIME
	
	ALPHA = 1.0 / Wi;
	
	BETA = Effective_Pol_Vis * ALPHA * ALPHA;

// 	std::cout << "ALPHA = " << ALPHA << endl;
// 	std::cout << "BETA = " << BETA << endl;
	
	// zeros out Grad_Vel
	for(int i = 0; i <= 1; i++)
	{
	  for(int k = 0; k <= 1; k++)
	  {
	    Grad_Vel(i,k) = 0.0;
	  }
	}
	
	// Computes Grad_Vel = El_Vel * Vel_Gdsf' (matlab's transpose notation)
	for(int i = 0; i <= 1; i++)
	{
	  for(int l = 0; l <= 1; l++)
	  {
	    for(int k = 0; k <= Vel_Npe-1; k++)
	    {
	      Grad_Vel(i,l) += El_Vel(i,k) * Vel_Gdsf(l,k); // Gdsf transpose
	    }
	  }
	}

	// The invere of gradu times its determinant
	Det_Inv_Grad_Vel(0,0) = Grad_Vel(1,1);
	Det_Inv_Grad_Vel(0,1) = -Grad_Vel(0,1);
	Det_Inv_Grad_Vel(1,0) = -Grad_Vel(1,0);
	Det_Inv_Grad_Vel(1,1) = Grad_Vel(0,0);
	
	Det_Grad_Vel = Grad_Vel(0,0) * Grad_Vel(1,1) - Grad_Vel(0,1) * Grad_Vel(1,0);
	
	Trace_Grad_Vel = Grad_Vel(0,0) + Grad_Vel(1,1);

	Coeff = 1 - TIMESTEP * Trace_Grad_Vel + TIMESTEP * TIMESTEP * Det_Grad_Vel;
	
	// (I - k * det(gradu)*(gradu Inv)) transpose
	TempMat1(0,0) = 1 - TIMESTEP * Det_Inv_Grad_Vel(0,0);
	TempMat1(1,0) = - TIMESTEP * Det_Inv_Grad_Vel(0,1); // Note transpose is here
	TempMat1(0,1) = - TIMESTEP * Det_Inv_Grad_Vel(1,0); // and here (1,0) to (0,1), etc.
	TempMat1(1,1) = 1 - TIMESTEP * Det_Inv_Grad_Vel(1,1);
	
	// I - k * det(gradu)*(gradu Inv)
	TempMat2(0,0) = 1 - TIMESTEP * Det_Inv_Grad_Vel(0,0);
	TempMat2(0,1) = - TIMESTEP * Det_Inv_Grad_Vel(0,1);
	TempMat2(1,0) = - TIMESTEP * Det_Inv_Grad_Vel(1,0);
	TempMat2(1,1) = 1 - TIMESTEP * Det_Inv_Grad_Vel(1,1);
	
	for(int i = 0; i <= 1; i++)
	{
	  for(int k = 0; k <= 1; k++)
	  {
	    TempMat3(i,k) = 0.0;
	    TempMat4(i,k) = 0.0;
	  }
	}
	
	for(int i = 0; i <= 1; i++)
	{
	  for(int l = 0; l <= 1; l++)
	  {
	    for(int k = 0; k <= 1; k++)
	    {
	      TempMat3(i,l) += TempMat2(i,k) * El_Dep_Str(k,l);
	    }
	  }
	}

	for(int i = 0; i <= 1; i++)
	{
	  for(int l = 0; l <= 1; l++)
	  {
	    for(int k = 0; k <= 1; k++)
	    {
	      TempMat4(i,l) += 1.0/(Coeff * Coeff) * TempMat3(i,k) * TempMat1(k,l);
	    }
	  }
	}

	Inod = Pre_Nod(Ne,j) - 1;

	Str_New(Inod, 0) = 1.0 / (1.0 + TIMESTEP * ALPHA) * TempMat4(0,0);// 1.0 / (1.0 + TIMESTEP * ALPHA) * (TempMat4(0,0) + TIMESTEP * BETA);
	Str_New(Inod, 1) = 1.0 / (1.0 + TIMESTEP * ALPHA) * TempMat4(0,1);
	Str_New(Inod, 2) = 1.0 / (1.0 + TIMESTEP * ALPHA) * TempMat4(1,1);//1.0 / (1.0 + TIMESTEP * ALPHA) * (TempMat4(1,1) + TIMESTEP * BETA);

// 	std::cout << "XX = " << Str_New(Inod, 0) << ", XY = " << Str_New(Inod, 1) << ", YY = " << Str_New(Inod, 2) << endl;
	
      } // for j
    }
      
    // minus 1 for C++ indexing Ele 2: Ne = 1 and is Ne associated with an odd numbered element
    else if((Ne % 2 == 0) && (Ne <= 2 * NX - 2))
    {
      
      for (int j = 1; j <= Pre_Npe - 1; j++)  // loop over nodes 2 and 3 for the first row
      {
	Vel_Inod = Vel_Nod(Ne,j) - 1;

	alpha1 = Vel_Elxy(1,0) * Vel_Elxy(2,1) - Vel_Elxy(2,0) * Vel_Elxy(1,1); // x2 * y3 - x3 * y2
	alpha2 = Vel_Elxy(2,0) * Vel_Elxy(0,1) - Vel_Elxy(0,0) * Vel_Elxy(2,1); // x3 * y1 - x1 * y3
	alpha3 = Vel_Elxy(0,0) * Vel_Elxy(1,1) - Vel_Elxy(1,0) * Vel_Elxy(0,1); // x1 * y2 - x2 * y1
	
	// In general: 2 * Area = x2 * y3 - x3 * y2 + x3 * y1 - x1 * y3 + x1 * y2 - x2 * y1
	Two_Area = alpha1 + alpha2 + alpha3;
	x = Vel_Glxy(Vel_Inod, 0);
	y = Vel_Glxy(Vel_Inod, 1);
	x1 = Pre_Elxy(0,0);
	x2 = Pre_Elxy(1,0);
	x3 = Pre_Elxy(2,0);
	y1 = Pre_Elxy(0,1);
	y2 = Pre_Elxy(1,1);
	y3 = Pre_Elxy(2,1);
	
	Xi  = 1.0 / Two_Area * ((x - x3) * (y2 - y3) - (y - y3) * (x2 - x3));
	Eta = 1.0 / Two_Area * ((x - x1) * (y3 - y1) + (y - y1) * (x1 - x3));

	Shape2d(Xi, 
		Eta, 
		Vel_Elxy, 
		Vel_Npe, 
		VEL_FLAG, 
		Sf, 
		Vel_Gdsf, 
		DetJ);

	Shape2d(Xi, 
		Eta, 
		Pre_Elxy, 
		Pre_Npe, 
		PRE_FLAG, 
		Pre_Sf, 
		Gdsf, 
		DetJ);

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

	Cur_Ele = Ne + 1;
	
	// Compute the stress at the departure foot
	for (int i = 0; i <= Pre_Npe - 1; i++) 	// get global coordinates of local nodes of element NE
	{
	  Inod = Pre_Nod(Cur_Ele-1,i) - 1;	// Global node number (minus one for C++ indexing) of local node.
	  Str_Elxy(i,0) = Pre_Glxy(Inod,0);	// x-coordinate
	  Str_Elxy(i,1) = Pre_Glxy(Inod,1);	// y-coordinate

	  El_Str(0,i) = Str_Old(Inod,0);
	  El_Str(1,i) = Str_Old(Inod,1);
	  El_Str(2,i) = El_Str(1,i); // The stress tensor is symmetric
	  El_Str(3,i) = Str_Old(Inod,2);

	}

	alpha1 = Str_Elxy(1,0) * Str_Elxy(2,1) - Str_Elxy(2,0) * Str_Elxy(1,1); // x2 * y3 - x3 * y2
	alpha2 = Str_Elxy(2,0) * Str_Elxy(0,1) - Str_Elxy(0,0) * Str_Elxy(2,1); // x3 * y1 - x1 * y3
	alpha3 = Str_Elxy(0,0) * Str_Elxy(1,1) - Str_Elxy(1,0) * Str_Elxy(0,1); // x1 * y2 - x2 * y1
	
	Two_Area = alpha1 + alpha2 + alpha3;
	
	Y_New_Xi_Eta(0) = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(2,0)) * (Str_Elxy(1,1) - Str_Elxy(2,1)) - 
			  (Depart_Foot(1) - Str_Elxy(2,1)) * (Str_Elxy(1,0) - Str_Elxy(2,0)));
	
	Y_New_Xi_Eta(1) = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(0,0)) * (Str_Elxy(2,1) - Str_Elxy(0,1)) + 
			  (Depart_Foot(1) - Str_Elxy(0,1)) * (Str_Elxy(0,0) - Str_Elxy(2,0)));
	
	//The point may have not left the element but it did move so recompute everything at the new position
	Xi = Y_New_Xi_Eta(0);
	Eta = Y_New_Xi_Eta(1);

	Shape2d(Xi, 
		Eta, 
		Str_Elxy, 
		Pre_Npe, 	// Both stress and pressure are linear
		PRE_FLAG, 
		Dep_Sf, 	// output
		Dep_Gdsf, // output
		DetJ);	// output

	for(int i = 0; i <= 2; i++)
	{
	  It_Str(i) = 0.0;
	}

	for(int k = 0; k <= Pre_Npe - 1; k++)
	{
	  It_Str(0) += El_Str(0,k) * Dep_Sf(k);
	  It_Str(1) += El_Str(1,k) * Dep_Sf(k);
	  It_Str(2) += El_Str(3,k) * Dep_Sf(k); // Since the Str is symmetric
	}
	
	Depart_Str(0) = It_Str(0);
	Depart_Str(1) = It_Str(1);
	Depart_Str(2) = It_Str(2);
	
	// This is not the element departure stress but the Stress at the particular departed Guass point
	El_Dep_Str(0,0) = Depart_Str(0);
	El_Dep_Str(0,1) = Depart_Str(1);
	El_Dep_Str(1,0) = El_Dep_Str(0,1); // The stress tensor is symmetric
	El_Dep_Str(1,1) = Depart_Str(2);

	// zeros out Grad_Vel
	for(int i = 0; i <= 1; i++)
	{
	  for(int k = 0; k <= 1; k++)
	  {
	    Grad_Vel(i,k) = 0.0;
	  }
	}
	
	// Computes Grad_Vel = El_Vel * Vel_Gdsf' (matlab's transpose notation)
	for(int i = 0; i <= 1; i++)
	{
	  for(int l = 0; l <= 1; l++)
	  {
	    for(int k = 0; k <= Vel_Npe-1; k++)
	    {
	      Grad_Vel(i,l) += El_Vel(i,k) * Vel_Gdsf(l,k); // Gdsf transpose
	    }
	  }
	}

	// The invere of gradu times its determinant
	Det_Inv_Grad_Vel(0,0) = Grad_Vel(1,1);
	Det_Inv_Grad_Vel(0,1) = -Grad_Vel(0,1);
	Det_Inv_Grad_Vel(1,0) = -Grad_Vel(1,0);
	Det_Inv_Grad_Vel(1,1) = Grad_Vel(0,0);
	
	Det_Grad_Vel = Grad_Vel(0,0) * Grad_Vel(1,1) - Grad_Vel(0,1) * Grad_Vel(1,0);
	
	Trace_Grad_Vel = Grad_Vel(0,0) + Grad_Vel(1,1);

	Coeff = 1 - TIMESTEP * Trace_Grad_Vel + TIMESTEP * TIMESTEP * Det_Grad_Vel;
	
	// (I - k * det(gradu)*(gradu Inv)) transpose
	TempMat1(0,0) = 1 - TIMESTEP * Det_Inv_Grad_Vel(0,0);
	TempMat1(1,0) = - TIMESTEP * Det_Inv_Grad_Vel(0,1); // Note transpose is here
	TempMat1(0,1) = - TIMESTEP * Det_Inv_Grad_Vel(1,0); // and here (1,0) to (0,1), etc.
	TempMat1(1,1) = 1 - TIMESTEP * Det_Inv_Grad_Vel(1,1);
	
	// I - k * det(gradu)*(gradu Inv)
	TempMat2(0,0) = 1 - TIMESTEP * Det_Inv_Grad_Vel(0,0);
	TempMat2(0,1) = - TIMESTEP * Det_Inv_Grad_Vel(0,1);
	TempMat2(1,0) = - TIMESTEP * Det_Inv_Grad_Vel(1,0);
	TempMat2(1,1) = 1 - TIMESTEP * Det_Inv_Grad_Vel(1,1);
	
	for(int i = 0; i <= 1; i++)
	{
	  for(int k = 0; k <= 1; k++)
	  {
	    TempMat3(i,k) = 0.0;
	    TempMat4(i,k) = 0.0;
	  }
	}
	
	for(int i = 0; i <= 1; i++)
	{
	  for(int l = 0; l <= 1; l++)
	  {
	    for(int k = 0; k <= 1; k++)
	    {
	      TempMat3(i,l) += TempMat2(i,k) * El_Dep_Str(k,l);
	    }
	  }
	}

	for(int i = 0; i <= 1; i++)
	{
	  for(int l = 0; l <= 1; l++)
	  {
	    for(int k = 0; k <= 1; k++)
	    {
	      TempMat4(i,l) += TempMat3(i,k) * TempMat1(k,l);
	    }
	  }
	}
	
	for(int i = 0; i <= 1; i++)
	{
	  for(int j = 0; j <= 1; j++)
	  {
	    TempMat4(i,j) = 1.0/(Coeff * Coeff) * TempMat4(i,j);
	  }
	}

	Inod = Pre_Nod(Ne,j) - 1;

	Str_New(Inod, 0) = 1.0 / (1.0 + TIMESTEP * ALPHA) * TempMat4(0,0);// 1.0 / (1.0 + TIMESTEP * ALPHA) * (TempMat4(0,0) + TIMESTEP * BETA);
	Str_New(Inod, 1) = 1.0 / (1.0 + TIMESTEP * ALPHA) * TempMat4(0,1);
	Str_New(Inod, 2) = 1.0 / (1.0 + TIMESTEP * ALPHA) * TempMat4(1,1);//1.0 / (1.0 + TIMESTEP * ALPHA) * (TempMat4(1,1) + TIMESTEP * BETA);

// 	std::cout << "XX = " << Str_New(Inod, 0) << ", XY = " << Str_New(Inod, 1) << ", YY = " << Str_New(Inod, 2) << endl;
	
      } // for j
    }
      
    // If Ne is the first even element in a row or an odd element greater than those in the first row
    // The additional of the 2 * NX is to shift the first even elements value larger than 2 * NX - 1
    // Since this is modular arithmetic this shift has no impact on larger values of Ne.
    // or Ne is an odd element larger than the those in the first row.
    else if(((Ne + 1) % (2 * NX) == 2) || (Ne % 2 == 0))
    {
      JJ = Pre_Npe - 1; // Only solve at the 3rd node of each element
    
      Vel_Inod = Vel_Nod(Ne,JJ) - 1;

      alpha1 = Vel_Elxy(1,0) * Vel_Elxy(2,1) - Vel_Elxy(2,0) * Vel_Elxy(1,1); // x2 * y3 - x3 * y2
      alpha2 = Vel_Elxy(2,0) * Vel_Elxy(0,1) - Vel_Elxy(0,0) * Vel_Elxy(2,1); // x3 * y1 - x1 * y3
      alpha3 = Vel_Elxy(0,0) * Vel_Elxy(1,1) - Vel_Elxy(1,0) * Vel_Elxy(0,1); // x1 * y2 - x2 * y1
      
      // In general: 2 * Area = x2 * y3 - x3 * y2 + x3 * y1 - x1 * y3 + x1 * y2 - x2 * y1
      Two_Area = alpha1 + alpha2 + alpha3;
      x = Vel_Glxy(Vel_Inod, 0);
      y = Vel_Glxy(Vel_Inod, 1);
      x1 = Pre_Elxy(0,0);
      x2 = Pre_Elxy(1,0);
      x3 = Pre_Elxy(2,0);
      y1 = Pre_Elxy(0,1);
      y2 = Pre_Elxy(1,1);
      y3 = Pre_Elxy(2,1);
      
      Xi  = 1.0 / Two_Area * ((x - x3) * (y2 - y3) - (y - y3) * (x2 - x3));
      Eta = 1.0 / Two_Area * ((x - x1) * (y3 - y1) + (y - y1) * (x1 - x3));

      Shape2d(Xi, 
	      Eta, 
	      Vel_Elxy, 
	      Vel_Npe, 
	      VEL_FLAG, 
	      Sf, 
	      Vel_Gdsf, 
	      DetJ);

      Shape2d(Xi, 
	      Eta, 
	      Pre_Elxy, 
	      Pre_Npe, 
	      PRE_FLAG, 
	      Pre_Sf, 
	      Gdsf, 
	      DetJ);

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

      Cur_Ele = Ne + 1;
      
      // Compute the stress at the departure foot
      for (int i = 0; i <= Pre_Npe - 1; i++) 	// get global coordinates of local nodes of element NE
      {
	Inod = Pre_Nod(Cur_Ele-1,i) - 1;	// Global node number (minus one for C++ indexing) of local node.
	Str_Elxy(i,0) = Pre_Glxy(Inod,0);	// x-coordinate
	Str_Elxy(i,1) = Pre_Glxy(Inod,1);	// y-coordinate

	El_Str(0,i) = Str_Old(Inod,0);
	El_Str(1,i) = Str_Old(Inod,1);
	El_Str(2,i) = El_Str(1,i); // The stress tensor is symmetric
	El_Str(3,i) = Str_Old(Inod,2);

      }

      alpha1 = Str_Elxy(1,0) * Str_Elxy(2,1) - Str_Elxy(2,0) * Str_Elxy(1,1); // x2 * y3 - x3 * y2
      alpha2 = Str_Elxy(2,0) * Str_Elxy(0,1) - Str_Elxy(0,0) * Str_Elxy(2,1); // x3 * y1 - x1 * y3
      alpha3 = Str_Elxy(0,0) * Str_Elxy(1,1) - Str_Elxy(1,0) * Str_Elxy(0,1); // x1 * y2 - x2 * y1
      
      Two_Area = alpha1 + alpha2 + alpha3;
      
      Y_New_Xi_Eta(0) = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(2,0)) * (Str_Elxy(1,1) - Str_Elxy(2,1)) - 
			(Depart_Foot(1) - Str_Elxy(2,1)) * (Str_Elxy(1,0) - Str_Elxy(2,0)));
      
      Y_New_Xi_Eta(1) = 1.0 / Two_Area * ((Depart_Foot(0) - Str_Elxy(0,0)) * (Str_Elxy(2,1) - Str_Elxy(0,1)) + 
			(Depart_Foot(1) - Str_Elxy(0,1)) * (Str_Elxy(0,0) - Str_Elxy(2,0)));
      
      //The point may have not left the element but it did move so recompute everything at the new position
      Xi = Y_New_Xi_Eta(0);
      Eta = Y_New_Xi_Eta(1);

      Shape2d(Xi, 
	      Eta, 
	      Str_Elxy, 
	      Pre_Npe, 	// Both stress and pressure are linear
	      PRE_FLAG, 
	      Dep_Sf, 	// output
	      Dep_Gdsf, // output
	      DetJ);	// output

      for(int i = 0; i <= 2; i++)
      {
	It_Str(i) = 0.0;
      }

      for(int k = 0; k <= Pre_Npe - 1; k++)
      {
	It_Str(0) += El_Str(0,k) * Dep_Sf(k);
	It_Str(1) += El_Str(1,k) * Dep_Sf(k);
	It_Str(2) += El_Str(3,k) * Dep_Sf(k); // Since the Str is symmetric
      }
      
      Depart_Str(0) = It_Str(0);
      Depart_Str(1) = It_Str(1);
      Depart_Str(2) = It_Str(2);
      
      // This is not the element departure stress but the Stress at the particular departed Guass point
      El_Dep_Str(0,0) = Depart_Str(0);
      El_Dep_Str(0,1) = Depart_Str(1);
      El_Dep_Str(1,0) = El_Dep_Str(0,1); // The stress tensor is symmetric
      El_Dep_Str(1,1) = Depart_Str(2);

      // zeros out Grad_Vel
      for(int i = 0; i <= 1; i++)
      {
	for(int k = 0; k <= 1; k++)
	{
	  Grad_Vel(i,k) = 0.0;
	}
      }
      
      // Computes Grad_Vel = El_Vel * Vel_Gdsf' (matlab's transpose notation)
      for(int i = 0; i <= 1; i++)
      {
	for(int l = 0; l <= 1; l++)
	{
	  for(int k = 0; k <= Vel_Npe-1; k++)
	  {
	    Grad_Vel(i,l) += El_Vel(i,k) * Vel_Gdsf(l,k); // Gdsf transpose
	  }
	}
      }

//       std::cout << "Vel_Gdsf = " << Vel_Gdsf << endl;
//       std::cout << "El_Vel = " << El_Vel << endl;
//       std::cout << "Gradu = " << Grad_Vel << endl;
      
      // The invere of gradu times its determinant
      Det_Inv_Grad_Vel(0,0) = Grad_Vel(1,1);
      Det_Inv_Grad_Vel(0,1) = -Grad_Vel(0,1);
      Det_Inv_Grad_Vel(1,0) = -Grad_Vel(1,0);
      Det_Inv_Grad_Vel(1,1) = Grad_Vel(0,0);
      
      Det_Grad_Vel = Grad_Vel(0,0) * Grad_Vel(1,1) - Grad_Vel(0,1) * Grad_Vel(1,0);
      
      Trace_Grad_Vel = Grad_Vel(0,0) + Grad_Vel(1,1);

      Coeff = 1 - TIMESTEP * Trace_Grad_Vel + TIMESTEP * TIMESTEP * Det_Grad_Vel;
      
      // (I - k * det(gradu)*(gradu Inv)) transpose
      TempMat1(0,0) = 1 - TIMESTEP * Det_Inv_Grad_Vel(0,0);
      TempMat1(1,0) = - TIMESTEP * Det_Inv_Grad_Vel(0,1); // Note transpose is here
      TempMat1(0,1) = - TIMESTEP * Det_Inv_Grad_Vel(1,0); // and here (1,0) to (0,1), etc.
      TempMat1(1,1) = 1 - TIMESTEP * Det_Inv_Grad_Vel(1,1);
      
      // I - k * det(gradu)*(gradu Inv)
      TempMat2(0,0) = 1 - TIMESTEP * Det_Inv_Grad_Vel(0,0);
      TempMat2(0,1) = - TIMESTEP * Det_Inv_Grad_Vel(0,1);
      TempMat2(1,0) = - TIMESTEP * Det_Inv_Grad_Vel(1,0);
      TempMat2(1,1) = 1 - TIMESTEP * Det_Inv_Grad_Vel(1,1);
      
//       std::cout << "F = " << TempMat2 << endl;
      
      for(int i = 0; i <= 1; i++)
      {
	for(int k = 0; k <= 1; k++)
	{
	  TempMat3(i,k) = 0.0;
	  TempMat4(i,k) = 0.0;
	}
      }
      
      for(int i = 0; i <= 1; i++)
      {
	for(int l = 0; l <= 1; l++)
	{
	  for(int k = 0; k <= 1; k++)
	  {
	    TempMat3(i,l) += TempMat2(i,k) * El_Dep_Str(k,l);
	  }
	}
      }

      for(int i = 0; i <= 1; i++)
      {
	for(int l = 0; l <= 1; l++)
	{
	  for(int k = 0; k <= 1; k++)
	  {
	    TempMat4(i,l) += TempMat3(i,k) * TempMat1(k,l);
	  }
	}
      }

      for(int i = 0; i <= 1; i++)
	{
	  for(int j = 0; j <= 1; j++)
	  {
	    TempMat4(i,j) = 1.0/(Coeff * Coeff) * TempMat4(i,j);
	  }
	}
      
      Inod = Pre_Nod(Ne, JJ) - 1;

      Str_New(Inod, 0) = 1.0 / (1.0 + TIMESTEP * ALPHA) * (TempMat4(0,0) + TIMESTEP * BETA);
      Str_New(Inod, 1) = 1.0 / (1.0 + TIMESTEP * ALPHA) * TempMat4(0,1);
      Str_New(Inod, 2) = 1.0 / (1.0 + TIMESTEP * ALPHA) * (TempMat4(1,1) + TIMESTEP * BETA);
      
//       std::cout << "XX = " << Str_New(Inod, 0) << ", XY = " << Str_New(Inod, 1) << ", YY = " << Str_New(Inod, 2) << endl;

    }// end if
  } // for NE
}