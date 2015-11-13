#include "InitialGuess.h"

void InitialGuess(int *ProcNodes, 
		  int Npe,
		  int Vel_Nnm,
		  int *ProcEles,
		  E_ISDM Vel_Nod,
		  E_ISDM Pre_Nod,
		  int TotalNumEles,
		  int myid,
		  double *PreviousSoln, 
		  E_V & x_VP)
{
  int GlobalNode, Error;
  
  double Values[1]; //, Temp_Values[2 * Vel_Nnm + Pre_Nnm];
  int Indices[1]; //, Temp_Indices[2 * Vel_Nnm + Pre_Nnm];

  for(int i = 0; i < TotalNumEles; i++)
  {
    if(ProcEles[i] == myid)
    {
      for(int j = 0; j < Npe; j++)
      {
	GlobalNode = Vel_Nod(i,j) - 1;
	if(ProcNodes[GlobalNode] == myid)
	{

	  // Horizontal Vel
	  Values[0] = PreviousSoln[GlobalNode];
	  Indices[0] = GlobalNode;
	  
	  Error = x_VP.ReplaceGlobalValues(1, Values, Indices);
	  
// 	  if(Error != 0)
// 	  {
// 	    std::cout << "Broke 1" << endl;
// 	    Error = 0;
// 	    
// 	    std::cout << "GlobalNode = " << GlobalNode << endl;
// 	    std::cout << "Element = " << i << endl;
// 	    std::cout << "Local Node = " << j << endl;
// 	    std::cout << "PreviousSoln[GlobalNode] = " << PreviousSoln[GlobalNode] << endl;
// 	    std::cout << "Values[0] = " << Values[0] << endl;
// 	    std::cout << "Indices[0] = " << Indices[0] << endl;
// 	    
// // 	    std::cin >> Error;
// 	  }
	  
	  // Vertical Vel
	  Values[0] = PreviousSoln[GlobalNode + Vel_Nnm];
	  Indices[0] = GlobalNode + Vel_Nnm;
	  
	  Error = x_VP.ReplaceGlobalValues(1, Values, Indices);
	  
// 	  if(Error != 0)
// 	  {
// 	    std::cout << "Broke 2" << endl;
// 	    Error = 0;
// 	    
// 	    std::cout << "GlobalNode = " << GlobalNode << endl;
// 	    std::cout << "Element = " << i << endl;
// 	    std::cout << "Local Node = " << j << endl;
// 	    std::cout << "PreviousSoln[GlobalNode] = " << PreviousSoln[GlobalNode] << endl;
// 	    std::cout << "Values[0] = " << Values[0] << endl;
// 	    std::cout << "Indices[0] = " << Indices[0] << endl;
// 	    
// // 	    std::cin >> Error;
// 	  }
	  
	  // Pressure
	  if((j == 0) | (j == 1) | (j == 2))
	  {
	    GlobalNode = Pre_Nod(i,j) - 1 + 2 * Vel_Nnm;
	    
	    Values[0] = PreviousSoln[GlobalNode];
	    Indices[0] = GlobalNode;
	    
	    Error = x_VP.ReplaceGlobalValues(1, Values, Indices);
	    
// 	    if(Error != 0)
// 	    {
// 	      std::cout << "Broke 3" << endl;
// 	      Error = 0;
// 	      
// 	      std::cout << "GlobalNode = " << GlobalNode << endl;
// 	      std::cout << "Element = " << i << endl;
// 	      std::cout << "Local Node = " << j << endl;
// 	      std::cout << "PreviousSoln[GlobalNode] = " << PreviousSoln[GlobalNode] << endl;
// 	      std::cout << "Values[0] = " << Values[0] << endl;
// 	      std::cout << "Indices[0] = " << Indices[0] << endl;
// 	      
// // 	      std::cin >> Error;
// 	    }
	  }
	}
      }
    }
  }
}