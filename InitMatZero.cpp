// Initializes the solution matrix to zero for every entry that will have a 
// value incerted.

#include "InitMatZero.h"

void InitMatZero(int Vel_Npe, 
		  int Pre_Npe,
		  int Nem,
		  int Vel_Nnm,
		  int Pre_Nnm,
		  E_ISDM Vel_Nod, 
		  E_ISDM Pre_Nod,
		  E_ISDM Vel_Nod_BC_Hor,
		  E_ISDM Vel_Nod_BC_Ver,
		  E_ISDM Pre_Nod_BC,
		  E_ISDV All_Proc_Nodes_VP,
		  E_ISDV My_Proc_Eles,
		  int myid,
		  E_CM & A) {

  // There are, at most, 6 horizontal vel, 6 vertical vel and 3 pre nodes that a processor 
  // could own.  As such there are at most 15 components computed per element andd at most 15
  // values "entered" in a row of the matrix per element.
  double *Values;
  Values = new double[15];
  int *Indices;
  Indices = new int[15];

  // Pre rows
  // There are at most 12 vel components per row
  double *P_Values = new double[12];
  int *P_Indices = new int[12];

  int ii, jj, kk, Error, IndexCounter;
  
  for (int Ne = 0; Ne < Nem; Ne++) 		// loop over all the elements
  {
    // If a processor has a node in the element then proceed
    if(My_Proc_Eles(Ne) == myid)
    {

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!! UPPER MATRICES START HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      // loop over nodes
      for (int i = 0; i <= Vel_Npe-1; i++)
      {
	
	// Initialize the index counter
	IndexCounter = 0;

	// The row of the global matrix
	ii = Vel_Nod(Ne,i) - 1;

	// If the node belongs to the processor then proceed
	if(All_Proc_Nodes_VP(ii) == myid)
	{
	  // global node ii is not a BC and hence the row is not an "identity" row
	  if(Vel_Nod_BC_Hor(Ne,i) != 1)
	  {
	    // Loop over the nodes in the element again for the column entries
	    for(int j = 0; j <= Vel_Npe-1; j++)
	    {
	      // The column of the global matrix
	      jj = Vel_Nod(Ne,j) - 1;

	      // Only record a zero if the column is not associated with a BC because
	      // if the column is associated with a BC then nothing will be enteried into
	      // the matrix but rather the RHS will get updated.
	      if(Vel_Nod_BC_Hor(Ne,j) != 1) // && Vel_Nod_BC_Hor(Ne,j) != 2)
	      {
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!! A11 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      // The horizontal velocity contribution to the horizontal velocity.
		Values[IndexCounter] = 0.0;
		
		Indices[IndexCounter] = jj;

		IndexCounter++;
		
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!! A12 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      // The vertical velocity contribution to the horizontal velocity
		Values[IndexCounter] = 0.0;
		
		Indices[IndexCounter] = jj + Vel_Nnm;
		
		IndexCounter++;
	      }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!! A12 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      // The vertical velocity contribution to the horizontal velocity
// 	      if(Vel_Nod_BC_Ver(Ne,j) != 1)
// 	      {
// 		Values[IndexCounter] = 0.0;
// 		
// 		Indices[IndexCounter] = jj + Vel_Nnm;
// 		
// 		IndexCounter++;
// 	      }
	    } // vel j

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!! B1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    // The pressure contribution to the horizontal velocity
	    for (int k = 0; k <= Pre_Npe-1; k++)
	    {
	      kk = Pre_Nod(Ne,k) - 1;

// 	      if(Pre_Nod_BC(Ne,k) != 1)
// 	      {

	      Values[IndexCounter] = 0.0;

	      Indices[IndexCounter] = kk + 2 * Vel_Nnm;
	      
	      IndexCounter++;
// 	      }
	    } // pre k

	    //Build Global Matrix A one row at a time
	    Error = A.InsertGlobalValues(ii, IndexCounter, Values, Indices);
	    
// 	    if(Error != 0 && myid == 1)
// 	    {
// 	      std::cout << "InitMat Zero Broke 1" << endl;
// 	      std::cout << "myid =" << myid << endl;
// 	      std::cout << "ii = " << ii << endl;
// 	      std::cout << "jj = " << jj << endl;
// 	      std::cout << "All_Proc_Nodes_VP(jj) = " << All_Proc_Nodes_VP(jj) << endl;
// 	      Error = 0;
// 	    }
	    
	  } // if ii != 1
	  else
	  {

	    // If the row is associated with a BC then it is an "identity" row
	    Values[0] = 1.0;

	    Indices[0] = ii;
	    
	    IndexCounter = 1;
	    
	    Error = A.InsertGlobalValues(ii, IndexCounter, Values, Indices);
	    
// 	    if(Error != 0)
// 	    {
// 	      std::cout << "InitMat Zero Broke 2" << endl;
// 	      std::cout << "myid =" << myid << endl;
// 	      std::cout << "ii = " << ii << endl;
// 	      std::cout << "jj = " << jj << endl;
// 	      std::cout << "All_Proc_Nodes_VP(ii) = " << All_Proc_Nodes_VP(ii) << endl;
// 	      Error = 0;
// 	    }
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
	
	// Everything that follows is logically the same as above but the following
	// conditional statement is fundamentally different.
	if(All_Proc_Nodes_VP(ii + Vel_Nnm) == myid)
	{
	  if(Vel_Nod_BC_Ver(Ne,i) != 1) // global node ii is not an EBC
	  {
	    // VEL COLUMNS
	    for(int j = 0; j <= Vel_Npe-1; j++)
	    {
	      // The column of the global matrix
	      jj = Vel_Nod(Ne,j) - 1;
	      
	      if(Vel_Nod_BC_Hor(Ne,j) != 1)// && Vel_Nod_BC_Hor(Ne,j) != 2)
	      {
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!! A21 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Values[IndexCounter] = 0.0;

		Indices[IndexCounter] = jj;
		
		IndexCounter++;
		
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!! A22 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Values[IndexCounter] = 0.0;

		Indices[IndexCounter] = jj + Vel_Nnm;
		
		IndexCounter++;
	      }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!! A22 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// 	      if(Vel_Nod_BC_Ver(Ne,j) != 1)
// 	      {
// 		Values[IndexCounter] = 0.0;
// 
// 		Indices[IndexCounter] = jj + Vel_Nnm;
// 		
// 		IndexCounter++;
// 	      }
	    } // vel j

	    // PRE COLUMNS CONT
	    for (int k = 0; k <= Pre_Npe-1; k++)
	    {
	      kk = Pre_Nod(Ne,k) - 1;

// 		if(Pre_Nod_BC(Ne,k) != 1)
// 		{

	      Values[IndexCounter] = 0.0;

	      Indices[IndexCounter] = kk + 2 * Vel_Nnm;

	      IndexCounter++;
// 		}
	    } // pre k

	    Error = A.InsertGlobalValues(ii + Vel_Nnm, IndexCounter, Values, Indices);
	    
// 	    if(Error != 0)
// 	    {
// 	      std::cout << "InitMat Zero Broke 3" << endl;
// 	      std::cout << "myid =" << myid << endl;
// 	      std::cout << "ii = " << ii << endl;
// 	      std::cout << "jj = " << jj << endl;
// 	      std::cout << "All_Proc_Nodes_VP(jj) = " << All_Proc_Nodes_VP(jj) << endl;
// 	      Error = 0;
// 	    }
	  } // if ii != 1
	  else // Essential BC on this row
	  {

	    Values[0] = 1.0;

	    Indices[0] = ii + Vel_Nnm;
	    
	    IndexCounter = 1;
	    
	    Error = A.InsertGlobalValues(ii + Vel_Nnm, IndexCounter, Values, Indices);
	    
// 	    if(Error != 0)
// 	    {
// 	      std::cout << "InitMat Zero Broke 4" << endl;
// 	      std::cout << "myid =" << myid << endl;
// 	      std::cout << "ii = " << ii << endl;
// 	      std::cout << "jj = " << jj << endl;
// 	      std::cout << "All_Proc_Nodes_VP(jj) = " << All_Proc_Nodes_VP(jj) << endl;
// 	      Error = 0;
// 	    }
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
	  // Pressure BCs are disabled because specifiying a pressure essential BC
	  // over-determines the system of equations that I am solving
// 	  if(Pre_Nod_BC(Ne,k) != 1)
// 	  {
	  // VEL COLUMNS
	  for(int j = 0; j <= Vel_Npe - 1; j++)
	  {
	    
	    jj = Vel_Nod(Ne,j) - 1;
	    
	    if(Vel_Nod_BC_Hor(Ne,j) != 1)// && Vel_Nod_BC_Hor(Ne,j) != 2)
	    {  
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!! B1 Trans !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    // The horizontal velocity contribution to the pressure
	      P_Values[IndexCounter] = 0.0;

	      P_Indices[IndexCounter] = jj;
	      
	      IndexCounter++;
	      
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!! B2 Trans !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    // The vertical velocity contribution to the pressure
	      P_Values[IndexCounter] = 0.0;

	      P_Indices[IndexCounter] = jj + Vel_Nnm;
	      
	      IndexCounter++;
	    }
	      
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!! B2 Trans !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    // The vertical velocity contribution to the pressure
// 	    if(Vel_Nod_BC_Ver(Ne,j) != 1)
// 	    {
// 	      P_Values[IndexCounter] = 0.0;
// 
// 	      P_Indices[IndexCounter] = jj + Vel_Nnm;
// 	      
// 	      IndexCounter++;
// 	    }
	  } // for j
	  
	  Error = A.InsertGlobalValues(kk + 2 * Vel_Nnm, IndexCounter, P_Values, P_Indices);
	  
// 	  if(Error != 0)
// 	  {
// 	    std::cout << "InitMat Zero Broke 5" << endl;
// 	    std::cout << "myid =" << myid << endl;
// 	    std::cout << "ii = " << ii << endl;
// 	    std::cout << "jj = " << jj << endl;
// 	    std::cout << "All_Proc_Nodes_VP(jj) = " << All_Proc_Nodes_VP(jj) << endl;
// 	    Error = 0;
// 	  }
// 	  }
	  // kk will never equal jj.  This BC condition only happens twice
// 	  else // if(Pre_Nod_BC(Ne,k) == 1 && All_Proc_Nodes_VP(kk + 2 * Vel_Nnm) == myid)
// 	  {
// 
// 	    P_Values[0] = 0.0;
// 	    
// 	    P_Indices[0] = kk + 2 * Vel_Nnm;
// 	    
// 	    IndexCounter = 1;
// 
// 	    Error = A.InsertGlobalValues(kk + 2 * Vel_Nnm, IndexCounter, P_Values, P_Indices);
// 	    
// 	    if(Error != 0)
// 	    {
// 	      std::cout << "InitMat Zero Broke 6" << endl;
// 	      std::cout << "myid =" << myid << endl;
// 	      std::cout << "ii = " << ii << endl;
// 	      std::cout << "jj = " << jj << endl;
// 	      std::cout << "All_Proc_Nodes_VP(jj) = " << All_Proc_Nodes_VP(jj) << endl;
// 	      Error = 0;
// 	    }
// 	  } // if kk
	}
      } // pre k
    } // If All_Proc_Eles == Ne
  } // for Ne
  
//   P_Values[0] = 0.0;
//   
//   IndexCounter = 1;
//   
//   for(int i = 0; i < Pre_Nnm; i++)
//   {
//     if(All_Proc_Nodes_VP(i + 2 * Vel_Nnm) == myid)
//     {
//       P_Indices[0] = i + 2 * Vel_Nnm;
// 
//       Error = A.InsertGlobalValues(i + 2 * Vel_Nnm, IndexCounter, P_Values, P_Indices);
//       
// //       if(Error != 0)
// //       {
// // 	std::cout << "InitMat Zero Broke 7" << endl;
// // 	Error = 0;
// //       }
//     }
//   }

  delete [] Values;
  delete [] Indices;
  delete [] P_Values;
  delete [] P_Indices;
  
}