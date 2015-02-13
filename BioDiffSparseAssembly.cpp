// Assembles the solution matrix and the RHS

// U_Feet is NNMX3 where each row is the global node number the departure foot is assocaited 
// with and the entires in the row are the x,y coordinates of the departure foot followed by 
// the ellemnt the foot is located in


#include "BioSparseAssembly.h"

void BioDiffSparseAssembly(double DIFF_COEFF,
			    double TIMESTEP,
			    int NX,
			    int NY,
			    int NGP,
			    int Vel_Npe,
			    int Nem, 
			    int N_TRI_QUAD,
			    int Vel_Nnm,
			    E_ISDM Vel_Nod, 
			    E_ISDM Bio_Nod_BC,
			    E_SDM Vel_Glxy,
			    E_SDM Tri_Quad_Pt, 
			    E_SDV Tri_Quad_Wt,
			    E_SDM GAUSPT,
			    E_SDM GAUSWT,
			    double *Bio_Old,
			    E_ISDV All_Proc_Nodes_Bio,
			    E_ISDV My_Proc_Eles,
			    int myid,
			    E_CM & A_Bio, 
			    E_V & b_Bio) {
  
  E_SDM Elxy, Gdsf;
  E_SDV Sf, Ele_Bio;
  double Essen_Bio;

  // Other variables
  int ii, jj, Inod, IndexCounter, Error;//, Gauss_Pt_Num, Cur_Ele;
  double Xi, Eta, DetJ, Const;
//   double x, y;
  double Bio_temp;

  //Vel
  // There are, at most, 6 horizontal vel, 6 vertical vel and 3 pre components per row
  double *Values = new double[15];
  int *Indices = new int[15];
  // int *Indices2 = new int[15];
  
//   Nat_Bio.Size(2);

  // RHS
  double b_Values[2]; //, Temp_Values[Vel_Nnm];
  int b_Indices[2]; //, Temp_Indices[Vel_Nnm];

  // Initialize
  // Velocity terms
  Elxy.Shape(Vel_Npe,2);
  Sf.Size(Vel_Npe);  		// value of shape functions at (Xi,Eta)
  Gdsf.Shape(2,Vel_Npe); 		// derivatives w.r.t. global cooridinates
  Ele_Bio.Size(Vel_Npe);

  DIFF_COEFF = 0.0001;
  
  for (int Ne = 0; Ne <= Nem - 1; Ne++) 		// loop over all the elements
  {
    if(My_Proc_Eles(Ne) == myid)
    {

      for (int i = 0; i <= Vel_Npe - 1; i++) 	// get global coordinates of local nodes of element Ne
      {
	Inod = Vel_Nod(Ne,i) - 1;	// Global node number (minus one for indeXing) of local node.
	Elxy(i,0) = Vel_Glxy(Inod,0);	// x-coordinate of the velcity
	Elxy(i,1) = Vel_Glxy(Inod,1);	// y-coordinate
	Ele_Bio(i) = Bio_Old[Inod];
      }

      for (int Ni = 0; Ni < N_TRI_QUAD; Ni++) // loop over quadrature points
      {
	Xi  = Tri_Quad_Pt(Ni,0);
	Eta = Tri_Quad_Pt(Ni,1);

	Shape2d(Xi, 
		Eta, 
		Elxy,
		Vel_Npe, 
		2, 		// vel_flag = 2
		Sf, 		// output
		Gdsf,		// output
		DetJ);		// output

	Bio_temp = 0.0;

	for(int k = 0; k < Vel_Npe; k++)
	{
	  Bio_temp += Sf(k) * Ele_Bio(k);
	}

	// Constant from change of variables, scaled and multiplied by the Gauss weight, see page 559
	Const = 0.5 * Tri_Quad_Wt(Ni) * DetJ;  // These two are the same and can be reduced to one term

	for (int i = 0; i <= Vel_Npe-1; i++)
	{
	  
	  // Initialize the index counter
	  IndexCounter = 0;

	  // The row of the global matrix
	  ii = Vel_Nod(Ne,i) - 1;
	  
	  if(All_Proc_Nodes_Bio(ii) == myid)
	  {
	    if(Bio_Nod_BC(Ne,i) != 1) // global node ii is not a BC
	    {

	      //RHS Vector
	      b_Values[0] = Const * Bio_temp * Sf(i);

	      b_Indices[0] = ii;

	      Error = b_Bio.SumIntoGlobalValues(1, b_Values, b_Indices);
	      
	      if(Error != 0)
	      {
		std::cout << "Broke Bio 2" << endl;
		Error = 0;
	      }

	      for(int j = 0; j <= Vel_Npe-1; j++)
	      {

		// The column of the global matrix
		jj = Vel_Nod(Ne,j) - 1;

		if(Bio_Nod_BC(Ne,j) != 1) // global node jj is not a BC  Doesn't matter if the column is a Natural BC
		{

  		Values[IndexCounter] = Const * (Sf(i) * Sf(j) + DIFF_COEFF * TIMESTEP * (Gdsf(0,i) * Gdsf(0,j) + Gdsf(1,i) * Gdsf(1,j)));

		Indices[IndexCounter] = jj;

		IndexCounter++;

		} // if jj != 1
		else // (Bio_Nod_BC(Ne,j) == 1) // column of non-EC row is an EC
		{

		  // RHS = RHS - A(ii,jj) * EBC
    // 	      x = Vel_Glxy(jj,0);
    // 	      y = Vel_Glxy(jj,1);
		  
    // 	      BioEssenBoundary2d(x,
    // 				y,
    // 				Essen_Bio);
		  Essen_Bio = 0.0;

		  b_Values[0] = 0.0; // -Const * Essen_Bio * (Sf(i) * Sf(j) + DIFF_COEFF * TIMESTEP * (Gdsf(0,i) * Gdsf(0,j) + Gdsf(1,i) * Gdsf(1,j)));

		  b_Indices[0] = ii;

		  Error = b_Bio.SumIntoGlobalValues(1, b_Values, b_Indices);
		  
		  if(Error != 0)
		  {
		    std::cout << "Broke Bio 4" << endl;
		    Error = 0;
		  }
		} // if jj = 1
	      } // vel j

	      //Build Global Matrix A one row at a time
	      Error = A_Bio.SumIntoGlobalValues(ii, IndexCounter, Values, Indices);
	      
	      if(Error != 0)
	      {
		std::cout << "Broke Bio 7" << endl;
		Error = 0;
	      }
	      
	    } // if ii != 1
	    else // if(Bio_Nod_BC(Ne,i) == 1) // Essential BC on this row
	    {
	      // RHS = EBC;
    // 	  x = Vel_Glxy(ii,0);
    // 	  y = Vel_Glxy(ii,1);

    // 	  BioEssenBoundary2d(x,
    // 			    y,
    // 			    Essen_Bio);
	      
	      Essen_Bio = 0.0; // / Shared_Nodes_Bio(ii);
	      
	      b_Values[0] = Essen_Bio;// EBCs is not defined yet.

	      b_Indices[0] = ii;

	      Error = b_Bio.ReplaceGlobalValues(1, b_Values, b_Indices);
	      
	      if(Error != 0)
	      {
		std::cout << "Broke Bio 8" << endl;
		Error = 0;
	      }

	      for(int j = 0; j <= Vel_Npe-1; j++)
	      {
		
		// The column of the global matrix
		jj = Vel_Nod(Ne,j) - 1;

		if(ii == jj)
		{

		  Values[0] = 1.0; // / Shared_Nodes_Bio(ii);

		  Indices[0] = ii;
		  
		  Error = A_Bio.ReplaceGlobalValues(ii, 1, Values, Indices);
		  
		  if(Error != 0)
		  {
		    std::cout << "Broke Bio 9" << endl;
		    Error = 0;
		  }
		}
	      } // vel j
	    } // if ii for ii = 1
	  }
	} // vel i
      } // for Ni
    } // All_Proc_Eles(Ne)
  } // for Ne

  delete [] Values;
  delete [] Indices;
  
}