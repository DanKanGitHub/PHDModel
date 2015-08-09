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
		    int Str_Npe,
		    int Nem, 
		    int N_TRI_QUAD,
		    int Vel_Nnm,
		    int Pre_Nnm,
		    int Str_Nnm,
		    E_ISDM Vel_Nod,
		    E_ISDM Pre_Nod,
		    E_ISDM Vel_Nod_BC_Hor,
		    E_ISDM Vel_Nod_BC_Ver,
		    E_ISDM Pre_Nod_BC, 
		    E_SDM Vel_Glxy, 
		    E_SDM Pre_Glxy,
		    int VEL_FLAG,
		    int PRE_FLAG, 
		    int STRESS_FLAG,
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
		    E_SDM GaussDepartFootx,
		    E_SDM GaussDepartFooty,
		    E_ISDM GaussDepartElement,
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
  E_SDV Dep_Sf;
  int ii, jj, kk, New_Ele, Inod, IndexCounter, Error;//, Gauss_Pt_Num, Cur_Ele;
  double Xi, Eta, Depart_a00, a00, a11, a22, DetJ, Const;
  double x, y;
  double Temp_A11_Global, Temp_A12_Global, Temp_A21_Global, Temp_A22_Global;
  double Temp_B1_Global, Temp_B2_Global; //, Temp_Shared;
  double Effective_Den, Effective_Vis, Rey_Num, Depart_Rey_Num; // Effective_Sol_Vis, Effective_Pol_Vis, 
  double P_Zero, U_Zero;
  
  //Vel
  // There are, at most, 6 horizontal vel, 6 vertical vel and 3 pre components per row
  double *Values1 = new double[15];
  double *Values2 = new double[15]; // We fill two rows at a time
  int *Indices = new int[15];

  Essen_Vel.Size(2);
  Nat_Vel.Size(2);
  
  // Pre rows
  // There are 12 vel components per row
  double *P_Values = new double[12]; // We fill two rows at a time
  int *P_Indices = new int[12];
  
  // RHS
  double b_Values[2]; //, Temp_Values[2 * Vel_Nnm + Pre_Nnm];
  int b_Indices[2]; //, Temp_Indices[2 * Vel_Nnm + Pre_Nnm];

  El_Vel.Shape(2,Vel_Npe);
  El_Bio.Size(Vel_Npe);
  It_Vel.Size(2);
  Depart_Vel.Size(2);

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
  El_Str.Shape(4,Str_Npe);
  Gp_Str.Size(Str_Npe);			// Stress at a Gauss point
  
  // Other terms
  Dep_Gdsf.Shape(2,Vel_Npe);
  Dep_Sf.Size(Vel_Npe);

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
      }
      
      if(STRESS_FLAG == 1)
      {
	for (int i = 0; i <= Str_Npe-1; i++) 	// get global coordinates of local nodes of element Ne
	{
	  Pre_Inod = Pre_Nod(Ne,i)-1;
	  El_Str(0,i) = Str(Pre_Inod,0);
	  El_Str(1,i) = Str(Pre_Inod,1);
	  El_Str(2,i) = Str(Pre_Inod,1); // The Stress is symmetric
	  El_Str(3,i) = Str(Pre_Inod,2);
	}
      }
      else
      {
	for (int i = 0; i <= Str_Npe-1; i++) 	// get global coordinates of local nodes of element Ne
	{
	  Vel_Inod = Vel_Nod(Ne,i) - 1;
	  El_Str(0,i) = Str(Vel_Inod,0);
	  El_Str(1,i) = Str(Vel_Inod,1);
	  El_Str(2,i) = Str(Vel_Inod,1); // The Stress is symmetric
	  El_Str(3,i) = Str(Vel_Inod,2);
	}
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
	if(STRESS_FLAG == 1)
	{
	  for(int j = 0; j <= Str_Npe - 1; j++)
	  {
	    Gp_Str(0) += El_Str(0,j) * Pre_Sf(j);
	    Gp_Str(1) += El_Str(1,j) * Pre_Sf(j);
	    Gp_Str(2) += El_Str(3,j) * Pre_Sf(j);
	  }
	}
	else
	{
	  for(int j = 0; j <= Str_Npe - 1; j++)
	  {
	    Gp_Str(0) += El_Str(0,j) * Vel_Sf(j);
	    Gp_Str(1) += El_Str(1,j) * Vel_Sf(j);
	    Gp_Str(2) += El_Str(3,j) * Vel_Sf(j);
	  }
	}

	New_Ele = GaussDepartElement(Ne, Ni);
	
// 	std::cout << "Here ok?" << endl;
	
	for (int k = 0; k <= Vel_Npe-1; k++) // get global coordinates of local nodes of element NE
	{
	  Inod = Vel_Nod(New_Ele - 1,k) - 1;		// Global node number of local node.
	  Vel_Elxy(k,0) = Vel_Glxy(Inod,0);	// x-coordinate of te new element
	  Vel_Elxy(k,1) = Vel_Glxy(Inod,1);	// y-coordinate of the new element
    
	  El_Vel(0,k) = Vel_Old[Inod];
	  El_Vel(1,k) = Vel_Old[Inod + Vel_Nnm];
	  
	  El_Bio(k) = Bio[Inod];
	}
	
	Shape2d(GaussDepartFootx(Ne, Ni), 
		GaussDepartFooty(Ne, Ni), 
		Vel_Elxy, 
		Vel_Npe, 
		VEL_FLAG, 
		Dep_Sf, 	// output
		Dep_Gdsf, 	// output
		DetJ);		// output

	// Zero out Iterative velocity
	for(int k = 0; k <= 1; k++)
	{
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
	
	Effective_Den = (1 - GaussPt_Bio) * Sol_Density + GaussPt_Bio * Poly_Density;
	
	P_Zero = Effective_Den * U_Zero * U_Zero;
	
	Effective_Vis = (Sol_Vis * (1 - GaussPt_Bio_Weight) + Poly_Vis * GaussPt_Bio_Weight) / (T_ZERO * P_Zero);

	Rey_Num = Effective_Den * U_Zero * L_ZERO / Effective_Vis;

	a00 = Rey_Num / TIMESTEP;
	a11 = RetardDivRelax * Effective_Vis; // Effective Non-Dimensional Solvent Viscosity
	a22 = a11 / 2.0;

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!! UPPER MATRICES START HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	for (int i = 0; i <= Vel_Npe-1; i++)
	{
	  
	  // Initialize the index counter
	  IndexCounter = 0;

	  // The row of the global matrix
	  ii = Vel_Nod(Ne,i) - 1;

	  if(All_Proc_Nodes_VP(ii) == myid)
	  {
	    if(Vel_Nod_BC_Hor(Ne,i) != 1) // global node ii is not a BC
	    {

	      //RHS Vector
	      b_Values[0] = Const * (Vel_Sf(i) * Depart_a00 * Depart_Vel(0) - Bio_Depart_Weight * (
		Vel_Gdsf(0,i) * Gp_Str(0) + Vel_Gdsf(1,i) * Gp_Str(1)));

	      b_Indices[0] = ii;

	      Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);

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

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!! A11 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		if(Vel_Nod_BC_Hor(Ne,j) != 1) // global node jj is not a BC  DOesn't matter if the column is a Natural BC
		{

		  Values1[IndexCounter] = Const * (a11 * Vel_Gdsf(0,i) * Vel_Gdsf(0,j)
		    + a22 * Vel_Gdsf(1,i) * Vel_Gdsf(1,j) + a00 * Vel_Sf(i) * Vel_Sf(j));
		  
		  Indices[IndexCounter] = jj;

		  IndexCounter++;

		} // if jj != 1

		
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

		if(Vel_Nod_BC_Ver(Ne,j) != 1) // global node jj is not a BC  DOesn't matter if the column is a Natural BC
		{

		  Values1[IndexCounter] = Const * a22 * Vel_Gdsf(0,i) * Vel_Gdsf(1,j);
		  
		  Indices[IndexCounter] = jj + Vel_Nnm;

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

		  b_Indices[0] = ii;
		  
		  Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		  if(Error != 0)
		  {
		    std::cout << "Broke 5" << endl;
		    Error = 0;
		  }
		} // if jj = 1
	      } // vel j

	      // PRE COLUMNS CONT
	      for (int k = 0; k <= Pre_Npe-1; k++)
	      {
		kk = Pre_Nod(Ne,k) - 1;

		if(Pre_Nod_BC(Ne,k) != 1)
		{
		  // From B1
		  Values1[IndexCounter] = -Const * Vel_Gdsf(0,i) * Pre_Sf(k);

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

		  Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		  
		  if(Error != 0)
		  {
		    std::cout << "Broke 6" << endl;
		    Error = 0;
		  }
		} // if kk
	      } // pre k

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

	      // RHS = EBC;
	      x = Vel_Glxy(ii,0);
	      y = Vel_Glxy(ii,1);

	      VelEssenBoundary2d(NX,
				  Vel_Nnm,
				  ii,
				  x,
				  y,
				  Essen_Vel);

	      // Top Row
	      b_Values[0] = Essen_Vel(0); // / Shared_Nodes_VP(ii); // * Temp_Shared;// EBCs is not defined yet.

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

		if(ii == jj)
		{

		  Values1[0] = 1.0; // / Shared_Nodes_VP(ii);

		  Indices[0] = ii;
		  
		  Error = A.ReplaceGlobalValues(ii, 1, Values1, Indices);
		  
		  if(Error != 0)
		  {
		    std::cout << "Broke 9" << endl;
		    Error = 0;
		  }
		}
	      } // vel j
	    } // if ii for ii = 1
	  }
	} // vel i

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!! LOWER MATRICES START HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

		if(Vel_Nod_BC_Hor(Ne,j) != 1)
		{

		  Values2[IndexCounter] = Const * a22 * Vel_Gdsf(0,j) * Vel_Gdsf(1,i);

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

		  b_Indices[0] = ii + Vel_Nnm;
		  
		  Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		  
		  if(Error != 0)
		  {
		    std::cout << "Broke 11" << endl;
		    Error = 0;
		  }
		} // if jj = 1

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!! A22 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		if(Vel_Nod_BC_Ver(Ne,j) != 1) // global node jj is not a BC  DOesn't matter if the column is a Natural BC
		{

		  Values2[IndexCounter] = Const * (a22 * Vel_Gdsf(0,i) * Vel_Gdsf(0,j)
		    + a11 * Vel_Gdsf(1,i) * Vel_Gdsf(1,j) + a00 * Vel_Sf(i) * Vel_Sf(j));

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

		  b_Indices[0] = ii + Vel_Nnm;
		  
		  Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		  
		  if(Error != 0)
		  {
		    std::cout << "Broke 12" << endl;
		    Error = 0;
		  }
		} // if jj = 1
	      } // vel j

	      // PRE COLUMNS CONT
	      for (int k = 0; k <= Pre_Npe-1; k++)
	      {
		kk = Pre_Nod(Ne,k) - 1;

		if(Pre_Nod_BC(Ne,k) != 1)
		{

		  Values2[IndexCounter] = -Const * Vel_Gdsf(1,i) * Pre_Sf(k);

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

		  Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		  
		  if(Error != 0)
		  {
		    std::cout << "Broke 13" << endl;
		    Error = 0;
		  }
		} // if kk
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

	      b_Values[0] = Essen_Vel(1); // / Shared_Nodes_VP(ii); // * Temp_Shared;// EBCs is not defined yet.

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
	    // VEL COLUMNS
	    for(int j = 0; j <= Vel_Npe - 1; j++)
	    {
	      
	      jj = Vel_Nod(Ne,j) - 1;
	      
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!! B1 Trans !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	      if(Vel_Nod_BC_Hor(Ne,j) != 1)
	      {
		// From B1
		P_Values[IndexCounter] = -Const * Vel_Gdsf(0,j) * Pre_Sf(k);

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

		b_Indices[0] = kk + 2 * Vel_Nnm;
		
		Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		
		if(Error != 0)
		{
		  std::cout << "Broke 17" << endl;
		  Error = 0;
		}
	      } // if jj

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!! B2 Trans !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	      if(Vel_Nod_BC_Ver(Ne,j) != 1)
	      {

		// From B2
		P_Values[IndexCounter] = -Const * Vel_Gdsf(1,j) * Pre_Sf(k);

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

		b_Indices[0] = kk + 2 * Vel_Nnm;
		
		Error = b.SumIntoGlobalValues(1, b_Values, b_Indices);
		
		if(Error != 0)
		{
		  std::cout << "Broke 18" << endl;
		  Error = 0;
		}
	      } // if jj
	    } // for j
	    Error = A.SumIntoGlobalValues(kk + 2 * Vel_Nnm, IndexCounter, P_Values, P_Indices);
	    
	    if(Error != 0)
	    {
	      std::cout << "Broke 19" << endl;
	      Error = 0;
	    }
	  }// if kk is mine
	} // pre k
      } // for Ni
    } // ele if
//     std::cout << "Element Number = " << Ne << endl;
  } // for Ne

  delete [] Values1;
  delete [] Values2;
  delete [] Indices;
  delete [] P_Values;
  delete [] P_Indices;
  
}