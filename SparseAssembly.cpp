// Assembles the solution matrix and the RHS

// U_Feet is NNMX3 where each row is the global node number the departure foot is assocaited 
// with and the entires in the row are the x,y coordinates of the departure foot followed by 
// the ellemnt the foot is located in

//   Effective_Sol_Vis = RETARD_TIME / RELAX_TIME * SOL_NEW_VIS / (T_ZERO * P_Zero);	//InitialStr and SparseAssembly
//   Effective_Pol_Vis = (1 - RETARD_TIME / RELAX_TIME) * POL_NEW_VIS / (T_ZERO * P_Zero); // Non-Dim correction is / (T_ZERO * P_Zero)
//   Effective_Vis = 0.2 * POL_NEW_VIS + 0.8 * SOL_NEW_VIS;	// Go into SparseAssembly
//   Rey_Num = U_Zero * L_ZERO / Effective_Vis * Effective_Den;	// SparseAssembly only	

#include "SparseAssembly.h"

void SparseAssembly(double Sol_Vis, 
		    double Poly_Vis, 
		    double Sol_Density,
		    double Poly_Density,
		    double T_ZERO,
		    double L_ZERO,
		    double TIMESTEP, 
		    int NX,
		    int NY,
		    int NGP,
		    int Vel_Npe, 
		    int Pre_Npe,
		    int Nem, 
		    int N_TRI_QUAD,
		    int Vel_Nnm,
		    int Pre_Nnm,
		    E_ISDM Vel_Nod,
		    E_ISDM Pre_Nod,
		    E_ISDM Vel_Nod_BC_Hor,
		    E_ISDM Vel_Nod_BC_Ver,
		    E_ISDM Pre_Nod_BC, 
		    E_SDM Vel_Glxy, 
		    E_SDM Pre_Glxy, 
		    E_ISDM Ele_Neigh,
		    int VEL_FLAG,
		    int PRE_FLAG, 
		    E_SDM Tri_Quad_Pt, 
		    E_SDV Tri_Quad_Wt,
		    E_SDM GAUSPT,
		    E_SDM GAUSWT,
		    E_SDM Str, 
		    double *Bio,
		    double *Vel_Old,
		    E_ISDV All_Proc_Nodes_VP,
		    E_ISDV My_Proc_Eles,
		    int myid,
		    E_CM & A, 
		    E_V & b) {
  
  // Velocity variables
  E_SDM Vel_Elxy, Vel_Gdsf;
  E_SDV Vel_Sf, Depart_Vel;
  E_SDM El_Vel;
  E_SDV It_Vel, Essen_Vel, Nat_Vel;
  int Vel_Inod;
  
  // Pressure variables
  E_SDM Pre_Elxy, Pre_Gdsf;
  E_SDV Pre_Sf;
  int Pre_Inod;
  double Essen_Pre;
  
  // Stress variables
  E_SDM El_Str; // Stress on a particular element
  E_SDV Gp_Str; // Stress at a particular Gauss point
  
  //Bio
  E_SDV El_Bio;
  double GaussPt_Bio, Depart_Bio, Bio_Depart_Weight, RetardDivRelax;
  double GaussPt_Bio_Weight;
  
  // Other variables
  E_SDM Dep_Gdsf;
  E_SDV Ini_Foot, Depart_Foot, Y_New_Xi_Eta, Dep_Sf;
  int ii, jj, kk, New_Ele, Inod, IndexCounter, Error;//, Gauss_Pt_Num, Cur_Ele;
  double Xi, Eta, Depart_a00, a00, a11, a22, DetJ, Const;
  double alpha1, alpha2, alpha3, Two_Area;
  double a1, a2, b1, b2, c1, c2, x_Foot, y_Foot, x, y;
  double Temp_A11_Global, Temp_A12_Global, Temp_A21_Global, Temp_A22_Global;
  double Temp_B1_Global, Temp_B2_Global; //, Temp_Shared;
  double Effective_Den, Effective_Vis, Rey_Num, Depart_Rey_Num; // Effective_Sol_Vis, Effective_Pol_Vis, 
  double P_Zero, U_Zero;
  
  //Vel
  // There are, at most, 6 horizontal vel, 6 vertical vel and 3 pre components per row
  double *Values1 = new double[15];
  double *Values2 = new double[15]; // We fill two rows at a time
  int *Indices = new int[15];
  // int *Indices2 = new int[15];
  
  Essen_Vel.Size(2);
  Nat_Vel.Size(2);
  
  // Pre rows
  // There are 12 vel components per row
  double *P_Values = new double[12]; // We fill two rows at a time
  int *P_Indices = new int[12];
  
  // RHS
  double b_Values[2]; //, Temp_Values[2 * Vel_Nnm + Pre_Nnm];
  int b_Indices[2]; //, Temp_Indices[2 * Vel_Nnm + Pre_Nnm];
  
  // Create data array
//   double Data[15];
//   double * DataPtr = Data;
  
  // Create Index array
//   int Index[15];
//   int * IndexPtr = Index;

  El_Vel.Shape(2,Vel_Npe);
  El_Bio.Size(Vel_Npe);
  It_Vel.Size(2);
  Ini_Foot.Size(2);
  Depart_Foot.Size(2);
  Depart_Vel.Size(2);
  Y_New_Xi_Eta.Size(2);

  // Initialize
  // Velocity terms
  Vel_Elxy.Shape(Vel_Npe,2);
  Vel_Sf.Size(Vel_Npe);  		// value of shape functions at (Xi,Eta)
  Vel_Gdsf.Shape(2,Vel_Npe); 		// derivatives w.r.t. global cooridinates
  
  // Pressure terms
  Pre_Elxy.Shape(Pre_Npe,2);
  Pre_Sf.Size(Pre_Npe);  		// value of shape functions at (Xi,Eta)
  Pre_Gdsf.Shape(2,Pre_Npe); 		// derivatives w.r.t. global cooridinates
  
  // Stress terms
  El_Str.Shape(4,3);
  Gp_Str.Size(3);			// Stress at a Gauss point
  
  // Other terms
  Dep_Gdsf.Shape(2,Vel_Npe);
  Dep_Sf.Size(Vel_Npe);

  // PDE coefffiecients
  // I recycle these variables in the calculation below
//   a00 = RE / TIMESTEP;
//   a11 = NEW_VIS;
//   a22 = a11 / 2.0;

  U_Zero = L_ZERO / T_ZERO;
  
  for (int Ne = 0; Ne <= Nem - 1; Ne++) 		// loop over all the elements
  {
    if(My_Proc_Eles(Ne) == myid)
    {
      for (int i = 0; i <= Vel_Npe - 1; i++) 	// get global coordinates of local nodes of element Ne
      {
	Vel_Inod = Vel_Nod(Ne,i) - 1;		// Global node number (minus one for indeXing) of local node.
	Vel_Elxy(i,0) = Vel_Glxy(Vel_Inod,0);	// x-coordinate of the velcity
	Vel_Elxy(i,1) = Vel_Glxy(Vel_Inod,1);	// y-coordinate
	
	El_Vel(0,i) = Vel_Old[Vel_Inod];
	El_Vel(1,i) = Vel_Old[Vel_Inod + Vel_Nnm];
	
	El_Bio(i) = Bio[Vel_Inod];
      }

      for (int i = 0; i <= Pre_Npe-1; i++) 	// get global coordinates of local nodes of element Ne
      {
	Pre_Inod = Pre_Nod(Ne,i)-1;		// Global node number (minus one for indeXing) of local node.
	Pre_Elxy(i,0) = Pre_Glxy(Pre_Inod,0);	// x-coordinate of the pressure
	Pre_Elxy(i,1) = Pre_Glxy(Pre_Inod,1);	// y-coordinate
	
	El_Str(0,i) = Str(Pre_Inod,0);
	El_Str(1,i) = Str(Pre_Inod,1);
	El_Str(2,i) = Str(Pre_Inod,1); // The Stress is symmetric
	El_Str(3,i) = Str(Pre_Inod,2);
      }

      for (int Ni = 0; Ni <= N_TRI_QUAD-1; Ni++) // loop over quadrature points
      {
	Xi  = Tri_Quad_Pt(Ni,0);
	Eta = Tri_Quad_Pt(Ni,1);

	Shape2d(Xi, 
		Eta, 
		Vel_Elxy, 
		Vel_Npe, 
		VEL_FLAG, 
		Vel_Sf, 		// output
		Vel_Gdsf,		// output
		DetJ);			// output

	Shape2d(Xi, 
		Eta, 
		Pre_Elxy, 
		Pre_Npe, 
		PRE_FLAG, 
		Pre_Sf, 		// output
		Pre_Gdsf,		// output
		DetJ);			// output
	  
	// Constant from change of variables, scaled and multiplied by the Gauss weight, see page 559
	Const = 0.5 * Tri_Quad_Wt(Ni) * DetJ;  // These two are the same and can be reduced to one term

	// Zero out data structure for the stress at the Gauss points
	for(int i = 0; i <= Pre_Npe - 1; i++)
	{
	  Gp_Str(i) = 0.0;
	}

	// Stress at a Guass point
	for(int j = 0; j <= Pre_Npe - 1; j++)
	{
	  Gp_Str(0) += El_Str(0,j) * Pre_Sf(j);
	  Gp_Str(1) += El_Str(1,j) * Pre_Sf(j);
	  Gp_Str(2) += El_Str(3,j) * Pre_Sf(j);
	}

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
	x_Foot = 1.0 / (b1 * c2 - c1 * b2) * (c2 * a1 - c1 * a2);
      
	y_Foot = 1.0 / (b1 * c2 - c1 * b2) * (b1 * a2 - b2 * a1);
	
	Ini_Foot(0) = x_Foot;
	Ini_Foot(1) = y_Foot;

	DepartureFoot(VEL_FLAG, 
		      Ne, 
		      Vel_Nnm,
		      TIMESTEP,
		      Vel_Old, 
		      Vel_Glxy, 
		      Vel_Npe, 
		      Vel_Nod, 
		      Ele_Neigh,
		      Ini_Foot,		// in x,y space
		      Depart_Foot,	// output in x,y space
		      New_Ele);		// output

	// Determine the element the Gauss point has moved into
	if(New_Ele-1 != Ne)
	{
	  for (int k = 0; k <= Vel_Npe-1; k++) // get global coordinates of local nodes of element NE
	  {
	    Inod = Vel_Nod(New_Ele-1,k) - 1;		// Global node number of local node.
	    Vel_Elxy(k,0) = Vel_Glxy(Inod,0);	// x-coordinate of te new element
	    Vel_Elxy(k,1) = Vel_Glxy(Inod,1);	// y-coordinate of the new element
      
	    El_Vel(0,k) = Vel_Old[Inod];
	    El_Vel(1,k) = Vel_Old[Inod + Vel_Nnm];
	    
	    El_Bio(k) = Bio[Inod];
	  }
	}

	// Convert the cartesian coordinates to barycentric coordinates
	alpha1 = Vel_Elxy(1,0) * Vel_Elxy(2,1) - Vel_Elxy(2,0) * Vel_Elxy(1,1); // x2 * y3 - x3 * y2
	alpha2 = Vel_Elxy(2,0) * Vel_Elxy(0,1) - Vel_Elxy(0,0) * Vel_Elxy(2,1); // x3 * y1 - x1 * y3
	alpha3 = Vel_Elxy(0,0) * Vel_Elxy(1,1) - Vel_Elxy(1,0) * Vel_Elxy(0,1); // x1 * y2 - x2 * y1
	
	Two_Area = alpha1 + alpha2 + alpha3;
	
	Y_New_Xi_Eta(0) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(2,0)) * (Vel_Elxy(1,1) - Vel_Elxy(2,1)) - 
			  (Depart_Foot(1) - Vel_Elxy(2,1)) * (Vel_Elxy(1,0) - Vel_Elxy(2,0)));
	
	Y_New_Xi_Eta(1) = 1.0 / Two_Area * ((Depart_Foot(0) - Vel_Elxy(0,0)) * (Vel_Elxy(2,1) - Vel_Elxy(0,1)) + 
			  (Depart_Foot(1) - Vel_Elxy(0,1)) * (Vel_Elxy(0,0) - Vel_Elxy(2,0)));

	//The point may have not left the element but it did move so recompute everything at the new position
	Xi = Y_New_Xi_Eta(0);
	Eta = Y_New_Xi_Eta(1);

	Shape2d(Xi, 
		Eta, 
		Vel_Elxy, 
		Vel_Npe, 
		VEL_FLAG, 
		Dep_Sf, 	// output
		Dep_Gdsf, 	// output
		DetJ);		// output

	// Zero out Iterative velocity
	for(int k = 0; k <= 1; k++)
	{
// 	  It_Vel(k) = 0.0;
	  Depart_Vel(k) = 0.0;
	}
	
	Depart_Bio = 0.0;
	
	// Compute iterative velocity
	for(int k = 0; k <= Vel_Npe - 1; k++)
	{
	  Depart_Vel(0) += El_Vel(0,k) * Dep_Sf(k);
	  Depart_Vel(1) += El_Vel(1,k) * Dep_Sf(k);
	  
	  Depart_Bio += El_Bio(k) * Dep_Sf(k);
	}
	
	// RHS at departure foot
	Bio_Depart_Weight = BioWeightFunc(Depart_Bio);

	RetardDivRelax = RetardationDividedByRelaxation(Bio_Depart_Weight);

	Effective_Den = (1 - Depart_Bio) * Sol_Density + Depart_Bio * Poly_Density;
	
	P_Zero = Effective_Den * U_Zero * U_Zero;
	
	Effective_Vis = (Sol_Vis * (1 - Bio_Depart_Weight) + Poly_Vis * Bio_Depart_Weight) / (T_ZERO * P_Zero);
	
	Depart_Rey_Num = Effective_Den * U_Zero * L_ZERO / Effective_Vis;
	
	Depart_a00 = Depart_Rey_Num / TIMESTEP;

	// LHS at Gauss Point
	GaussPt_Bio = 0.0;
	
	for(int k = 0; k <= Vel_Npe - 1; k++)
	{
	  GaussPt_Bio += El_Bio(k) * Vel_Sf(k);
	}
	
	GaussPt_Bio_Weight = BioWeightFunc(GaussPt_Bio);
	
	RetardDivRelax = RetardationDividedByRelaxation(GaussPt_Bio);
	
	Effective_Vis = (Sol_Vis * (1 - GaussPt_Bio_Weight) + Poly_Vis * GaussPt_Bio_Weight) / (T_ZERO * P_Zero);
	
	Effective_Den = (1 - GaussPt_Bio) * Sol_Density + GaussPt_Bio * Poly_Density;
	
	P_Zero = Effective_Den * U_Zero * U_Zero;
	
	Rey_Num = Effective_Den * U_Zero * L_ZERO / Effective_Vis;

// 	Effective_Den = 1.0; // 0.2 * DENSITY_BIO + 0.8 * DENSITY_WATER;
// 	P_Zero = Effective_Den * U_Zero * U_Zero;
// 	Effective_Vis = (Poly_Vis + Sol_Vis) / (T_ZERO * P_Zero) / 2;
// // 	Effective_Sol_Vis = 1079.0 / 1080.0 * Effective_Vis;
// // 	Effective_Pol_Vis = (1 - 1079.0 / 1080.0) * Effective_Vis;
// 	Rey_Num = Effective_Den * U_Zero * L_ZERO / Effective_Vis;
	
	a00 = Rey_Num / TIMESTEP;
	a11 = RetardDivRelax * Effective_Vis; // Effective Non-Dimensional Solvent Viscosity
	a22 = a11 / 2.0;
	
// 	Depart_a00 = a00; // Just testing
	
// 	Depart_Vel(0) = It_Vel(0);
// 	Depart_Vel(1) = It_Vel(1);

// 	if(myid == 0)
// 	{
// 	  std::cout << "Here 1" << endl;
// 	  int QWERTY;
// 	  std::cin >> QWERTY;
// 	}
	
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!! UPPER MATRICES START HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	for (int i = 0; i <= Vel_Npe-1; i++)
	{
	  
	  // Initialize the index counter
	  IndexCounter = 0;

	  // The row of the global matrix
	  ii = Vel_Nod(Ne,i) - 1;
	  
// 	  if(myid == 0)
// 	  {
// 	    
// 	    std::cout << "ii = " << ii << endl;
// 	    
// 	  }
	  
	  if(All_Proc_Nodes_VP(ii) == myid)
	  {
	    if(Vel_Nod_BC_Hor(Ne,i) != 1) // global node ii is not a BC
	    {
	      
// 	      if(myid == 0)
// 	      {
// 		std::cout << "Here 1c" << endl;
// 		int QWERTY;
// 		std::cin >> QWERTY;
// 	      }
	      
	      //RHS Vector
	      b_Values[0] = Const * (Vel_Sf(i) * Depart_a00 * Depart_Vel(0) - Bio_Depart_Weight * (
		Vel_Gdsf(0,i) * Gp_Str(0) + Vel_Gdsf(1,i) * Gp_Str(1)));

// 	      b_Values[0] = Const * (Vel_Sf(i) * Depart_a00 * Depart_Vel(0) - Vel_Gdsf(0,i) * Gp_Str(0) - Vel_Gdsf(1,i) * Gp_Str(1));

    // 	  b_Dense(ii) += b_Values[0];

	      b_Indices[0] = ii;

	      Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);

// 	      if(ii == 9 & myid == 1)
// 	      {
// 		std::cout << "ii = 6" << endl;
// 		std::cout << "Ne = " << Ne << endl;
// 		std::cout << "i = " << i << endl;
// 		std::cout << "b_Values[0] = " << b_Values[0] << endl;
// 		std::cout << "Const = " << Const << endl;
// 		std::cout << "Vel_Sf(i) = " << Vel_Sf(i) << endl;
// 		std::cout << "a00 = " << a00 << endl;
// 		std::cout << "Depart_Vel(0) = " << Depart_Vel(0) << endl;
// 		std::cout << "Vel_Gdsf(0,i) = " << Vel_Gdsf(0,i) << endl;
// 		std::cout << "Gp_Str(0) = " << Gp_Str(0) << endl;
// 		std::cout << "Vel_Gdsf(1,i) = " << Vel_Gdsf(1,i) << endl;
// 		std::cout << "Gp_Str(1) = " << Gp_Str(1) << endl;
// 		
// 		int GHLKJ;
// 		std::cin >> GHLKJ;
// 	      }
	      
	      if(Error != 0)
	      {
		std::cout << "Broke 2" << endl;
		Error = 0;
	      }

	      if(Vel_Nod_BC_Hor(Ne,i) == 2) // This row is a natural BC
	      {
		// Modify RHS only
		NatBoundary2d(GAUSPT, 
			      GAUSWT, 
			      NX,
			      NGP,
			      Vel_Glxy, 
			      Vel_Nod,
			      Ne,
			      Nat_Vel);	// output

		b_Values[0] = Nat_Vel(0);

    // 	    b_Dense(ii) += b_Values[0];

		b_Indices[0] = ii;

		Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		
		if(Error != 0)
		{
		  std::cout << "Broke 3" << endl;
		  Error = 0;
		}
	      }

	      for(int j = 0; j <= Vel_Npe-1; j++)
	      {

		// The column of the global matrix
		jj = Vel_Nod(Ne,j) - 1;
		
// 		if(myid == 0)
// 		{
// 		  
// 		  std::cout << "jj = " << jj << endl;
// 		  
// 		}
		
// 		if(myid == 0)
// 		{
// 		  std::cout << "Here 1b" << endl;
// 		  int QWERTY;
// 		  std::cin >> QWERTY;
// 		}

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!! A11 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
  // 	      if(All_Proc_Nodes_VP(jj) == myid)
  // 	      {
		
		  if(Vel_Nod_BC_Hor(Ne,j) != 1) // global node jj is not a BC  DOesn't matter if the column is a Natural BC
		  {

		    Values1[IndexCounter] = Const * (a11 * Vel_Gdsf(0,i) * Vel_Gdsf(0,j)
		      + a22 * Vel_Gdsf(1,i) * Vel_Gdsf(1,j) + a00 * Vel_Sf(i) * Vel_Sf(j));
		    
		    Indices[IndexCounter] = jj;
		    
      // 	      A_Dense(ii,jj) += Values1[IndexCounter];

		    IndexCounter++;

		  } // if jj != 1
  // 	      }
		
		if(Vel_Nod_BC_Hor(Ne,j) == 1) // column of non-EC row is an EC
		{

		  // RHS = RHS - A(ii,jj) * EBC
		  x = Vel_Glxy(jj,0);
		  y = Vel_Glxy(jj,1);
		  
		  VelEssenBoundary2d(NX,
				    Vel_Nnm,
				    jj,
				    x,
				    y,
				    Essen_Vel);

		  Temp_A11_Global = Const * (a11 * Vel_Gdsf(0,i) * Vel_Gdsf(0,j)
		      + a22 * Vel_Gdsf(1,i) * Vel_Gdsf(1,j) + a00 * Vel_Sf(i) * Vel_Sf(j));

		  b_Values[0] = -Temp_A11_Global * Essen_Vel(0);

    // 	      b_Dense(ii) += b_Values[0];

		  b_Indices[0] = ii;

		  Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		  if(Error != 0)
		  {
		    std::cout << "Broke 4" << endl;
		    Error = 0;
		  }
		} // if jj = 1

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!! A12 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
// 		if(myid == 0)
// 		{
// 		  std::cout << "Here 1a" << endl;
// 		  int QWERTY;
// 		  std::cin >> QWERTY;
// 		}
    
  // 	      if(All_Proc_Nodes_VP(jj + Vel_Nnm) == myid)
  // 	      {
		  if(Vel_Nod_BC_Ver(Ne,j) != 1) // global node jj is not a BC  DOesn't matter if the column is a Natural BC
		  {

		    Values1[IndexCounter] = Const * a22 * Vel_Gdsf(0,i) * Vel_Gdsf(1,j);
		    
		    Indices[IndexCounter] = jj + Vel_Nnm;

      // 	      A_Dense(ii,jj + Vel_Nnm) += Values1[IndexCounter];

		    IndexCounter++;
		  } // if jj != 1
		  else // if(Vel_Nod_BC_Ver(Ne,j) == 1) // column of non-EC row is an EC
		  {

		    // RHS = RHS - A(ii,jj) * EBC
		    x = Vel_Glxy(jj,0);
		    y = Vel_Glxy(jj,1);
		    
		    VelEssenBoundary2d(NX,
				      Vel_Nnm,
				      jj,
				      x,
				      y,
				      Essen_Vel);

		    Temp_A12_Global = Const * a22 * Vel_Gdsf(0,i) * Vel_Gdsf(1,j);
		    
		    b_Values[0] = -Temp_A12_Global * Essen_Vel(1);

      // 	      b_Dense(ii) += b_Values[0];

		    b_Indices[0] = ii;
		    
		    Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		    if(Error != 0)
		    {
		      std::cout << "Broke 5" << endl;
		      Error = 0;
		    }
		  } // if jj = 1
  // 	      }
	      } // vel j

	      // PRE COLUMNS CONT
	      for (int k = 0; k <= Pre_Npe-1; k++)
	      {
		kk = Pre_Nod(Ne,k) - 1;
		
  // 	      if(All_Proc_Nodes_VP(kk + 2 * Vel_Nnm) == myid)
  // 	      {

		  if(Pre_Nod_BC(Ne,k) != 1)
		  {
		    // From B1
		    Values1[IndexCounter] = -Const * Vel_Gdsf(0,i) * Pre_Sf(k);

      // 	      A_Dense(ii, kk + 2 * Vel_Nnm) += Values1[IndexCounter];

		    Indices[IndexCounter] = kk + 2 * Vel_Nnm;
		    
		    IndexCounter++;
		  }
		  else // if(Pre_Nod_BC(Ne,k) == 1)
		  {

		    // RHS = RHS - A(ii,jj) * EBC
		    x = Pre_Glxy(kk,0);
		    y = Pre_Glxy(kk,1);
		    Essen_Pre = PreEssenBoundary2d(x,
						    y);
		    
		    // Top Row
		    Temp_B1_Global = -Const * Vel_Gdsf(0,i) * Pre_Sf(k);

		    b_Values[0] = -Temp_B1_Global * Essen_Pre;

		    b_Indices[0] = ii;
		    
      // 	      b_Dense(ii) += b_Values[0];
		    
		    Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		    
		    if(Error != 0)
		    {
		      std::cout << "Broke 6" << endl;
		      Error = 0;
		    }
		  } // if kk
  // 	      }
	      } // pre k
	      
    // 	  std::cout << "ii = " << ii << endl;
    // 	  for(int i = 0; i < IndexCounter; i++)
    // 	  {
    // 	    
    // 	   std::cout << "A(" << ii << "," << Indices[i] << ") = " << Values1[i] << endl; 
    // 	    
    // 	  }
    // 	  
    // 	  int GHF;
    // 	  std::cin >> GHF;
	      
	      //Build Global Matrix A one row at a time
	      Error = A.SumIntoGlobalValues(ii, IndexCounter, Values1, Indices);
	      
	      if(Error != 0)
	      {
		std::cout << "Broke 7" << endl;
		Error = 0;
	      }
	      
	    } // if ii != 1
	    else // if(Vel_Nod_BC_Hor(Ne,i) == 1) // Essential BC on this row
	    {
	      
// 	      if(myid == 0)
// 	      {
// 		std::cout << "myid = " << myid << endl;
// 		std::cout << "Here 1d" << endl;
// 		std::cout << "ii = " << ii << endl;
// 	      }
	      
	      // RHS = EBC;
	      x = Vel_Glxy(ii,0);
	      y = Vel_Glxy(ii,1);

	      VelEssenBoundary2d(NX,
				  Vel_Nnm,
				  ii,
				  x,
				  y,
				  Essen_Vel);
	      
  // 	    if(Shared_Nodes_VP(ii) != 1.0)
  // 	    {
  // 	      
  // // 	      Temp_Shared = 1.0 / Shared_Nodes_VP(ii);
  // 	      std::cout << " Shared_Nodes_VP = " << Shared_Nodes_VP(ii) << endl;
  // 	      std::cout << "ii = " << ii << endl;
  // 
  // 	    }
  // 	    else
  // 	    {
  // 	      
  // 	      Temp_Shared = 1.0;
  // 	      
  // 	    }
	      
	      // Top Row
	      b_Values[0] = Essen_Vel(0); // / Shared_Nodes_VP(ii); // * Temp_Shared;// EBCs is not defined yet.

    // 	  b_Dense(ii) = b_Values[0];

	      b_Indices[0] = ii;

	      Error = b.ReplaceGlobalValues(1, b_Values, b_Indices);
	      
	      if(Error != 0)
	      {
		std::cout << "Broke 8" << endl;
		Error = 0;
	      }

	      //The following is now done in InitMatZero.cpp
	      for(int j = 0; j <= Vel_Npe-1; j++)
	      {
		
		// The column of the global matrix
		jj = Vel_Nod(Ne,j) - 1;

  // 	      if(All_Proc_Nodes_VP(jj) == myid)
  // 	      {
		
		  if(ii == jj)
		  {

		    Values1[0] = 1.0; // / Shared_Nodes_VP(ii);
      // 	      A_Dense(ii, jj) = 1.0;
		    Indices[0] = ii;
		    
		    Error = A.ReplaceGlobalValues(ii, 1, Values1, Indices);
		    
		    if(Error != 0)
		    {
		      std::cout << "Broke 9" << endl;
		      Error = 0;
		    }
		  }
  // 	      }
	      } // vel j
	    } // if ii for ii = 1
	  }
	} // vel i

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!! LOWER MATRICES START HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
// 	if(myid == 0)
// 	{
// 	  std::cout << "Here 2" << endl;
// 	  
// // 	  int QWERTY;
// // 	  std::cin >> QWERTY;
// 	}
	
	for (int i = 0; i <= Vel_Npe-1; i++)
	{
	  
	  // Initialize the index counter
	  IndexCounter = 0;

	  // The row of the global matrix
	  ii = Vel_Nod(Ne,i) - 1;
	  
	  if(All_Proc_Nodes_VP(ii + Vel_Nnm) == myid)
	  {
	    if(Vel_Nod_BC_Ver(Ne,i) != 1) // global node ii is not a BC
	    {
	      // The Bio_Depart_Weight scales the stress contribution by the amount of bio-film in the element
	      b_Values[0] = Const * (Vel_Sf(i) * Depart_a00 * Depart_Vel(1) - Bio_Depart_Weight * (
		Vel_Gdsf(0,i) * Gp_Str(1) + Vel_Gdsf(1,i) * Gp_Str(2)));
	      
// 	      b_Values[0] = Const * (Vel_Sf(i) * Depart_a00 * Depart_Vel(1) - Vel_Gdsf(0,i) * Gp_Str(1) - Vel_Gdsf(1,i) * Gp_Str(2));

    // 	  b_Dense(ii + Vel_Nnm) += b_Values[0];

	      b_Indices[0] = ii + Vel_Nnm;

	      Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
	      
	      if(Error != 0)
	      {
		std::cout << "Broke 10" << endl;
		Error = 0;
	      }

	      // VEL COLUMNS
	      for(int j = 0; j <= Vel_Npe-1; j++)
	      {

		// The column of the global matrix
		jj = Vel_Nod(Ne,j) - 1;
	      
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!! A21 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
  // 	      if(All_Proc_Nodes_VP(jj) == myid)
  // 	      {
		
		  if(Vel_Nod_BC_Hor(Ne,j) != 1)
		  {

		    Values2[IndexCounter] = Const * a22 * Vel_Gdsf(0,j) * Vel_Gdsf(1,i);
		    
      // 	      A_Dense(ii + Vel_Nnm,jj) += Values2[IndexCounter];
		    
		    Indices[IndexCounter] = jj;
		    
		    IndexCounter++;

		  } // if jj != 1
		  else // if(Vel_Nod_BC_Hor(Ne,j) == 1) // column of non-EC row is an EC
		  {

		    x = Vel_Glxy(jj,0);
		    y = Vel_Glxy(jj,1);
		    
		    VelEssenBoundary2d(NX,
				      Vel_Nnm,
				      jj,
				      x,
				      y,
				      Essen_Vel);

		    Temp_A21_Global = Const * a22 * Vel_Gdsf(0,j) * Vel_Gdsf(1,i);

		    b_Values[0] = -Temp_A21_Global * Essen_Vel(0);
		    
      // 	      b_Dense(ii + Vel_Nnm) += b_Values[0];

		    b_Indices[0] = ii + Vel_Nnm;
		    
		    Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		    
		    if(Error != 0)
		    {
		      std::cout << "Broke 11" << endl;
		      Error = 0;
		    }
		  } // if jj = 1
  // 	      }
		
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!! A22 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  // 	      if(All_Proc_Nodes_VP(jj + Vel_Nnm) == myid)
  // 	      {
		
		  if(Vel_Nod_BC_Ver(Ne,j) != 1) // global node jj is not a BC  DOesn't matter if the column is a Natural BC
		  {

		    Values2[IndexCounter] = Const * (a22 * Vel_Gdsf(0,i) * Vel_Gdsf(0,j)
		      + a11 * Vel_Gdsf(1,i) * Vel_Gdsf(1,j) + a00 * Vel_Sf(i) * Vel_Sf(j));
		    
      // 	      A_Dense(ii + Vel_Nnm,jj + Vel_Nnm) += Values2[IndexCounter];
		    
		    Indices[IndexCounter] = jj + Vel_Nnm;
		    
		    IndexCounter++;
		  } // if jj != 1
		  else // if(Vel_Nod_BC_Ver(Ne,j) == 1) // column of non-EC row is an EC
		  {

		    x = Vel_Glxy(jj,0);
		    y = Vel_Glxy(jj,1);
		    
		    VelEssenBoundary2d(NX,
				      Vel_Nnm,
				      jj,
				      x,
				      y,
				      Essen_Vel);

		    Temp_A22_Global = Const * (a22 * Vel_Gdsf(0,i) * Vel_Gdsf(0,j)
			+ a11 * Vel_Gdsf(1,i) * Vel_Gdsf(1,j) + a00 * Vel_Sf(i) * Vel_Sf(j));
		    
		    b_Values[0] = -Temp_A22_Global * Essen_Vel(1);
		    
      // 	      b_Dense(ii + Vel_Nnm) += b_Values[0];

		    b_Indices[0] = ii + Vel_Nnm;
		    
		    Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		    
		    if(Error != 0)
		    {
		      std::cout << "Broke 12" << endl;
		      Error = 0;
		    }
		  } // if jj = 1
  // 	      }
	      } // vel j

	      // PRE COLUMNS CONT
	      for (int k = 0; k <= Pre_Npe-1; k++)
	      {
		kk = Pre_Nod(Ne,k) - 1;

  // 	      if(All_Proc_Nodes_VP(kk + 2 * Vel_Nnm) == myid)
  // 	      {
		
		  if(Pre_Nod_BC(Ne,k) != 1)
		  {

		    Values2[IndexCounter] = -Const * Vel_Gdsf(1,i) * Pre_Sf(k);
		    
      // 	      A_Dense(ii + Vel_Nnm, kk + 2 * Vel_Nnm) += Values2[IndexCounter];
		    
		    Indices[IndexCounter] = kk + 2 * Vel_Nnm;

		    IndexCounter++;
		  }
		  // This BC condition only happens twice; first for Ele 1 node 1 and second for Ele 2 node 1.
		  else // if(Pre_Nod_BC(Ne,k) == 1)
		  {

		    x = Pre_Glxy(kk,0);
		    y = Pre_Glxy(kk,1);
		    Essen_Pre = PreEssenBoundary2d(x,
						    y);

		    Temp_B2_Global = -Const * Vel_Gdsf(1,i) * Pre_Sf(k);
		    
		    b_Values[0] = -Temp_B2_Global * Essen_Pre;

		    b_Indices[0] = ii + Vel_Nnm;
		    
      // 	      b_Dense(ii + Vel_Nnm) += b_Values[0];
		    
		    Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		    
		    if(Error != 0)
		    {
		      std::cout << "Broke 13" << endl;
		      Error = 0;
		    }
		  } // if kk
  // 	      }
	      } // pre k
	      
	      //Build Global Matrix A one row at a time
	      Error = A.SumIntoGlobalValues(ii + Vel_Nnm, IndexCounter, Values2, Indices);
	      
	      if(Error != 0)
	      {
		std::cout << "Broke 14" << endl;
		Error = 0;
	      }
	    } // if ii != 1
	    else // if(Vel_Nod_BC_Ver(Ne,i) == 1) // Essential BC on this row
	    {
	      // RHS = EBC;
	      x = Vel_Glxy(ii,0);
	      y = Vel_Glxy(ii,1);

	      VelEssenBoundary2d(NX,
				  Vel_Nnm,
				  ii,
				  x,
				  y,
				  Essen_Vel);
	      
  // 	    if(Shared_Nodes_VP(ii) != 1.0)
  // 	    {
  // 	      
  // 	      Temp_Shared = 1.0 / Shared_Nodes_VP(ii);
  // 	      
  // 	    }
  // 	    else
  // 	    {
  // 	      
  // 	      Temp_Shared = 1.0;
  // 	      
  // 	    }

	      b_Values[0] = Essen_Vel(1); // / Shared_Nodes_VP(ii); // * Temp_Shared;// EBCs is not defined yet.

    // 	  b_Dense(ii + Vel_Nnm) = b_Values[0];

	      b_Indices[0] = ii + Vel_Nnm;
	      
	      Error = b.ReplaceGlobalValues(1, b_Values, b_Indices);
	      
	      if(Error != 0)
	      {
		std::cout << "Broke 15" << endl;
		Error = 0;
	      }

	      for(int j = 0; j <= Vel_Npe-1; j++)
	      {
		
		// The column of the global matrix
		jj = Vel_Nod(Ne,j) - 1;

		if(ii == jj)
		{

		  Values2[0] = 1.0; // / Shared_Nodes_VP(ii);
		  
    // 	      A_Dense(ii + Vel_Nnm, jj + Vel_Nnm) = 1.0;
		  
		  Indices[0] = ii + Vel_Nnm;
		  
		  Error = A.ReplaceGlobalValues(ii + Vel_Nnm, 1, Values2, Indices);
		  
		  if(Error != 0)
		  {
		    std::cout << "Broke 16" << endl;
		    Error = 0;
		  }
		}
	      } // vel j
	    } // if ii for ii = 1
	  }
	} // vel i, part 2
	
	// Now do B1 transpose and B2 transpose
	for (int k = 0; k <= Pre_Npe - 1; k++)
	{
	  kk = Pre_Nod(Ne,k) - 1;

	  IndexCounter = 0;
	  
	  if(All_Proc_Nodes_VP(kk + 2 * Vel_Nnm) == myid)
	  {
	    // Pressure BC disabled
// 	    if(Pre_Nod_BC(Ne,k) != 1)
// 	    {

	    // VEL COLUMNS
	    for(int j = 0; j <= Vel_Npe - 1; j++)
	    {
	      
	      jj = Vel_Nod(Ne,j) - 1;
	      
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!! B1 Trans !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      
// 	      if(All_Proc_Nodes_VP(jj) == myid)
// 	      {
	      
		if(Vel_Nod_BC_Hor(Ne,j) != 1)
		{
		  // From B1
		  P_Values[IndexCounter] = -Const * Vel_Gdsf(0,j) * Pre_Sf(k);

    // 	      A_Dense(kk + 2 * Vel_Nnm, jj) += P_Values[IndexCounter];
		  
		  P_Indices[IndexCounter] = jj;
		  
		  IndexCounter++;

		}
		else // if(Vel_Nod_BC_Hor(Ne,j) == 1)
		{

		  // RHS = RHS - A(i,j) * EC
		  x = Vel_Glxy(jj,0);
		  y = Vel_Glxy(jj,1);
		  
		  VelEssenBoundary2d(NX,
				    Vel_Nnm,
				    jj,
				    x,
				    y,
				    Essen_Vel);
		  
		  //For horizontal velocity contribution
		  Temp_B1_Global = -Const * Vel_Gdsf(0,j) * Pre_Sf(k);

		  b_Values[0] = -Temp_B1_Global * Essen_Vel(0);

    // 	      b_Dense(kk + 2 * Vel_Nnm) += b_Values[0];

		  b_Indices[0] = kk + 2 * Vel_Nnm;
		  
		  Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		  
		  if(Error != 0)
		  {
		    std::cout << "Broke 17" << endl;
		    Error = 0;
		  }
		} // if jj
// 	      }
	      
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!! B2 Trans !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// 	      if(All_Proc_Nodes_VP(jj + Vel_Nnm) == myid)
// 	      {
		if(Vel_Nod_BC_Ver(Ne,j) != 1)
		{

		  // From B2
		  P_Values[IndexCounter] = -Const * Vel_Gdsf(1,j) * Pre_Sf(k);

    // 	      A_Dense(kk + 2 * Vel_Nnm, jj + Vel_Nnm) += P_Values[IndexCounter];
		  
		  P_Indices[IndexCounter] = jj + Vel_Nnm;
		  
		  IndexCounter++;
		}
		else // if(Vel_Nod_BC_Ver(Ne,j) == 1)
		{

		  // RHS = RHS - A(i,j) * EC
		  x = Vel_Glxy(jj,0);
		  y = Vel_Glxy(jj,1);
		  
		  VelEssenBoundary2d(NX,
				    Vel_Nnm,
				    jj,
				    x,
				    y,
				    Essen_Vel);

		  //For vertical velocity contribution
		  Temp_B2_Global = -Const * Vel_Gdsf(1,j) * Pre_Sf(k);

		  b_Values[0] = -Temp_B2_Global * Essen_Vel(1);

    // 	      b_Dense(kk + 2 * Vel_Nnm) += b_Values[0];

		  b_Indices[0] = kk + 2 * Vel_Nnm;
		  
		  Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		  
		  if(Error != 0)
		  {
		    std::cout << "Broke 18" << endl;
		    Error = 0;
		  }
		} // if jj
// 	      }
	    } // for j
	    Error = A.SumIntoGlobalValues(kk + 2 * Vel_Nnm, IndexCounter, P_Values, P_Indices);
	    
	    if(Error != 0)
	    {
	      std::cout << "Broke 19" << endl;
	      Error = 0;
	    }
// 	    }
	    // kk will never equal jj.  This BC condition only happens twice
// 	    else // if(Pre_Nod_BC(Ne,k) == 1)
// 	    {
// 
// 	      x = Pre_Glxy(kk,0);
// 	      y = Pre_Glxy(kk,1);
// 	      Essen_Pre = PreEssenBoundary2d(x,
// 					      y);
// 
// 	      b_Values[0] = Essen_Pre;
// 
//     // 	  b_Dense(kk + 2 * Vel_Nnm) = b_Values[0];
// 
//     // 	  A_Dense(kk + 2 * Vel_Nnm, kk + 2 * Vel_Nnm) = 1.0;
// 	      
// 	      b_Indices[0] = kk + 2 * Vel_Nnm;
// 	      
// 	      Error = b.ReplaceGlobalValues(1, b_Values, b_Indices);
// 	      
// 	      if(Error != 0)
// 	      {
// 		std::cout << "Broke 20" << endl;
// 		std::cout << "myid =" << myid << endl;
// 		std::cout << "ii = " << ii << endl;
// 		std::cout << "kk = " << kk << endl;
// 		std::cout << "All_Proc_Nodes_VP(ii) = " << All_Proc_Nodes_VP(kk + 2 * Vel_Nnm) << endl;
// 		Error = 0;
// 	      }
// 
// 	      P_Values[0] = 1.0;
// 	      P_Indices[0] = kk + 2 * Vel_Nnm;
// 
// 	      Error = A.ReplaceGlobalValues(kk + 2 * Vel_Nnm, 1, P_Values, P_Indices);
// 	      
// 	      if(Error != 0)
// 	      {
// 		std::cout << "Broke 21" << endl;
// 		std::cout << "myid =" << myid << endl;
// 		std::cout << "kk = " << kk << endl;
// 		std::cout << "All_Proc_Nodes_VP(ii) = " << All_Proc_Nodes_VP(kk + 2 * Vel_Nnm) << endl;
// 		Error = 0;
// 	      }
// 	    } // if BC row
	  }// if kk is mine
	} // pre k
      } // for Ni
    } // ele if
  } // for Ne

//   P_Values[0] = 0.000001;
//   
//   IndexCounter = 1;
//   
//   for(int i = 1; i < Pre_Nnm; i++)// Starts at 1 because i = 0 is taken care of above
//   {
// 	
//     P_Indices[0] = i + 2 * Vel_Nnm;
// 
//     Error = A.ReplaceGlobalValues(i + 2 * Vel_Nnm, IndexCounter, P_Values, P_Indices);
//     
//     if(Error != 0)
//     {
//       std::cout << "Broke 22" << endl; 
//     }
//   }
  
  delete [] Values1;
  delete [] Values2;
  delete [] Indices;
  delete [] P_Values;
  delete [] P_Indices;
  
}