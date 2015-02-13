// Processor node Partition

#include "ProcNodePartitionClass.h"

ProcNodePartitionClass::ProcNodePartitionClass()
{

}

ProcNodePartitionClass::~ProcNodePartitionClass()
{

}

// My implementation of Prathish's method
// void ProcNodePartitionClass::ProcNodePartition(int myid,
// 					       int NumProc, 
// 					       int Vel_Nnm, 
// 					       int Pre_Nnm, 
// 					       int Vel_Npe,
// 					       int Pre_Npe,
// 					       E_ISDM Vel_Nod,
// 					       E_ISDM Pre_Nod,
// 					       int Nem)
// {
// 
//   int Total_Num_Nodes, Num_Nodes, Mod_Nodes, count;
// 
//   Total_Num_Nodes = 2 * Vel_Nnm + Pre_Nnm;
//   
//   // Partition the pressure nodes first.
//   Mod_Nodes = Pre_Nnm % NumProc;
//   
//   Num_Nodes = Pre_Nnm / NumProc; // Since NumProc starts at zero.
//   
//   int *Nodes_Per_Proc = new int[NumProc];
//   int *Pre_Node_Part = new int[Pre_Nnm];
//   
//   All_Proc_Nodes_VP.Size(Total_Num_Nodes);
//   All_Proc_Nodes_Bio.Size(Vel_Nnm);
//   
//   // Pressure Num of nodes partitiioned
//   for(int i = 0; i < NumProc; i++)
//   {
//     if(i < Mod_Nodes)
//     {
//       Nodes_Per_Proc[i] = Num_Nodes + 1;
//     }
//     else
//     {
//       Nodes_Per_Proc[i] = Num_Nodes;
//     }
//   }
//   
//   // Assign a check value  Might be unnecessary
// //   for(int i = 0; i < Total_Num_Nodes; i++)
// //   {
// //     All_Proc_Nodes_VP(i) = NumProc;
// //   }
// //   
// //   for(int i = 0; i < Vel_Nnm; i++)
// //   {
// //     All_Proc_Nodes_Bio(i) = NumProc;
// //   }
//   
//   count = 0;
//   
//   // Assign pressure nodes myids
//   for(int i = 0; i < NumProc; i++)
//   {
//     for(int j = 0; j < Nodes_Per_Proc[i]; j++)
//     {
//     
//       Pre_Node_Part[count] = i;
//       All_Proc_Nodes_VP(count + 2 * Vel_Nnm) = i;
//       count++;
//       
//     }
//   }
//   
//   // Find associated velocity nodes for previously assigned pressure nodes
//   
//   for(int i = 0; i < NumProc; i++)
//   {
//     for(int e = 0; e < Nem; e++)
//     {
//       for(int n = 0; n < Pre_Npe; n++)
//       {
// 	if(Pre_Node_Part[Pre_Nod(e,n) - 1] == i)
// 	{
// 	  
// 	  All_Proc_Nodes_VP(Vel_Nod(e,n) - 1) = i;
// 	  All_Proc_Nodes_VP(Vel_Nod(e,n) - 1 + Vel_Nnm) = i;
// 	  
// 	  All_Proc_Nodes_Bio(Vel_Nod(e,n) - 1) = i;
// 	  
// 	}
//       }
//     }
//   }
// 
//   int Node1, Node2, Node3, Node4, Node5, Node6, myidN1, myidN2, myidN3;
//   
//   // Assign non-vertex nodes
//   for(int i = 0; i < Nem; i++)
//   {
//     
//     Node1 = Vel_Nod(i,0) - 1;
//     Node2 = Vel_Nod(i,1) - 1;
//     Node3 = Vel_Nod(i,2) - 1;
//     Node4 = Vel_Nod(i,3) - 1;
//     Node5 = Vel_Nod(i,4) - 1;
//     Node6 = Vel_Nod(i,5) - 1;
//     
//     myidN1 = All_Proc_Nodes_VP(Node1);
//     myidN2 = All_Proc_Nodes_VP(Node2);
//     myidN3 = All_Proc_Nodes_VP(Node3);
//     
//     if(myidN1 <= myidN2)
//     {
//       
//       All_Proc_Nodes_VP(Node4) = myidN1;
//       All_Proc_Nodes_VP(Node4 + Vel_Nnm) = myidN1;
//       
//       All_Proc_Nodes_Bio(Node4) = myidN1;
// 
//     }
//     
//     if(myidN1 > myidN2)
//     {
//       
//       All_Proc_Nodes_VP(Node4) = myidN2;
//       All_Proc_Nodes_VP(Node4 + Vel_Nnm) = myidN2;
//       
//       All_Proc_Nodes_Bio(Node4) = myidN2;
// 
//     }
//     
//     if(myidN1 <= myidN3)
//     {
//       
//       All_Proc_Nodes_VP(Node6) = myidN1;
//       All_Proc_Nodes_VP(Node6 + Vel_Nnm) = myidN1;
//       
//       All_Proc_Nodes_Bio(Node6) = myidN1;
// 
//     }
//     
//     if(myidN1 > myidN3)
//     {
//       
//       All_Proc_Nodes_VP(Node6) = myidN3;
//       All_Proc_Nodes_VP(Node6 + Vel_Nnm) = myidN3;
//       
//       All_Proc_Nodes_Bio(Node6) = myidN3;
// 
//     }
//     
//     if(myidN2 <= myidN3)
//     {
//       
//       All_Proc_Nodes_VP(Node5) = myidN2;
//       All_Proc_Nodes_VP(Node5 + Vel_Nnm) = myidN2;
//       
//       All_Proc_Nodes_Bio(Node5) = myidN2;
// 
//     }
//     
//     if(myidN2 > myidN3)
//     {
//       
//       All_Proc_Nodes_VP(Node5) = myidN3;
//       All_Proc_Nodes_VP(Node5 + Vel_Nnm) = myidN3;
//       
//       All_Proc_Nodes_Bio(Node5) = myidN3;
// 
//     }
//   }
//   
//   Nodes_Per_Proc_VP = 0;
//   Nodes_Per_Proc_Bio = 0;
//   
//   // Count how many nodes "I" own
//   for(int i = 0; i < Total_Num_Nodes; i++)
//   {
//     if(All_Proc_Nodes_VP(i) == myid)
//     {
//       
//       Nodes_Per_Proc_VP++;
//       
//     }
//   }
//   
//   for(int i = 0; i < Vel_Nnm; i++)
//   {
//     if(All_Proc_Nodes_Bio(i) == myid)
//     {
//       
//       Nodes_Per_Proc_Bio++;
//       
//     }
//   }
//   
//   // Build local global node arrays
//   My_Proc_Nodes_VP.Size(Nodes_Per_Proc_VP);
//   My_Proc_Nodes_Bio.Size(Nodes_Per_Proc_Bio);
//   
//   int new_my_count, new_my_count_Bio;
//   
//   new_my_count = 0;
//   new_my_count_Bio = 0;
//   
//   //Assign local nodes to arrays
//   for(int i = 0; i < Total_Num_Nodes; i++)
//   {
//     if(All_Proc_Nodes_VP(i) == myid)
//     {
//       
//       My_Proc_Nodes_VP(new_my_count) = i;
//       new_my_count++;
//       
//     }
//   }
//   
//   for(int i = 0; i < Vel_Nnm; i++)
//   {
//     if(All_Proc_Nodes_Bio(i) == myid)
//     {
//       
//       My_Proc_Nodes_Bio(new_my_count_Bio) = i;
//       new_my_count_Bio++;
//       
//     }
//   }
//   
//   delete [] Nodes_Per_Proc;
//   delete [] Pre_Node_Part;
// }

// My method
void ProcNodePartitionClass::ProcNodePartition(int myid,
					       int NumProc, 
					       int Vel_Nnm, 
					       int Pre_Nnm, 
					       int Vel_Npe,
					       int Pre_Npe,
					       E_ISDM Vel_Nod,
					       E_ISDM Pre_Nod,
					       int Nem)
{

  int Total_Num_Nodes, Num_Nodes, Mod_Nodes;

  Total_Num_Nodes = 2 * Vel_Nnm + Pre_Nnm;
  
  // Partition the pressure nodes first.
  Mod_Nodes = Pre_Nnm % NumProc;
  
  Num_Nodes = Pre_Nnm / NumProc; // Since NumProc starts at zero.
  
  int *Nodes_Per_Proc = new int[NumProc];
  
  All_Proc_Nodes_VP.Size(Total_Num_Nodes);
  All_Proc_Nodes_Bio.Size(Vel_Nnm);
  
  for(int i = 0; i < NumProc; i++)
  {
    if(i < Mod_Nodes)
    {
      Nodes_Per_Proc[i] = Num_Nodes + 1;
    }
    else
    {
      Nodes_Per_Proc[i] = Num_Nodes;
    }
  }

  for(int i = 0; i < Total_Num_Nodes; i++)
  {
    All_Proc_Nodes_VP(i) = NumProc;
  }
  
  for(int i = 0; i < Vel_Nnm; i++)
  {
    All_Proc_Nodes_Bio(i) = NumProc;
  }
  
//   std::cout << "Here 1" << endl;
//   std::cout << "All_Proc_Nodes_VP = " << All_Proc_Nodes_VP << endl;
//   for(int i = 0; i < NumProc; i++)
//   {
//     std::cout << "Nodes_Per_Proc = " << Nodes_Per_Proc[i] << endl;
//   }
  
//   int QWERTY;
//   std::cin >> QWERTY;
  
  int EleNum = 0;
  int count, NpeNum, Global_Node_VP, Global_Node_Pre;
  NpeNum = 0;
  
  for(int i = 0; i < NumProc; i++)
  {
    
    count = 0;
    
//     if(myid = 1)
//     {
//       std::cout << "i = " << i << endl;
//     }
    
    while(count < Nodes_Per_Proc[i]) // count is indexed with start at zero and Nodes_Per_Proc starts counting at 1
    {
      
      if(NpeNum == Pre_Npe)
      {
	
	NpeNum = 0;
	EleNum++;
	
      }
      
//       std::cout << "Vel_Nod = " << Vel_Nod << endl;
//       std::cout << "Pre_Nod = " << Pre_Nod << endl;

      Global_Node_VP = Vel_Nod(EleNum,NpeNum) - 1;
      Global_Node_Pre = Pre_Nod(EleNum,NpeNum) - 1;

//       std::cout << "Global_Node_VP = " << Global_Node_VP << endl;
//       std::cout << "Global_Node_Pre = " << Global_Node_Pre << endl;
      
//       if(myid == 1)
//       {
// 	std::cout << "EleNum = " << EleNum << endl;
// 	std::cout << "NpeNum = " << NpeNum << endl;
//       }
      
//       std::cout << "All_Proc_Nodes_VP = " << All_Proc_Nodes_VP << endl;
      
//       int QWERTY;
//       std::cin >> QWERTY;

      if(All_Proc_Nodes_VP(Global_Node_VP) == NumProc)
      {
	All_Proc_Nodes_VP(Global_Node_Pre + 2 * Vel_Nnm) = i;	// Pressure
	All_Proc_Nodes_VP(Global_Node_VP) = i;			// Hor Vel
	All_Proc_Nodes_VP(Global_Node_VP + Vel_Nnm) = i;	// Vert Vel
	
	All_Proc_Nodes_Bio(Global_Node_VP) = i;
	
// 	if(myid == 1)
// 	{
// 	  std::cout << "count = " << count << endl;
	  
// 	  int QWERTY;
// 	  std::cin >> QWERTY;
// 	}
	
	count++; // Only count when a node is assigned
      }
      
      NpeNum++;
      
    } // while
  } // for

//   std::cout << "All_Proc_Nodes_VP = " << All_Proc_Nodes_VP << endl;
//   
//   int QWERTY;
//   std::cin >> QWERTY;
  
  int Node1, Node2, Node3, Node4, Node5, Node6, myidN1, myidN2, myidN3;
  int new_my_count, new_my_count_Bio;
  
  new_my_count = 0;
  new_my_count_Bio = 0;
  
  My_Proc_Eles.Size(Nem);
  
  for(int i = 0; i < Nem; i++)
  {
    
    My_Proc_Eles(i) = -1;
    
  }
  
//   std::cout << "My_Proc_Eles = " << My_Proc_Eles << endl;
  
  for(int i = 0; i < Nem; i++)
  {
    
    Node1 = Vel_Nod(i,0) - 1;
    Node2 = Vel_Nod(i,1) - 1;
    Node3 = Vel_Nod(i,2) - 1;
    Node4 = Vel_Nod(i,3) - 1;
    Node5 = Vel_Nod(i,4) - 1;
    Node6 = Vel_Nod(i,5) - 1;
    
    myidN1 = All_Proc_Nodes_VP(Node1);
    myidN2 = All_Proc_Nodes_VP(Node2);
    myidN3 = All_Proc_Nodes_VP(Node3);
    
    if(myid == myidN1 || myid == myidN2 || myid == myidN3)
    {
      
      My_Proc_Eles(i) = myid;
      
    }

//     if(i == 7)
//     {
//       std::cout << "Node1 = " << Node1 << endl;
//       std::cout << "Node2 = " << Node2 << endl;
//       std::cout << "Node3 = " << Node3 << endl;
//       
//       std::cout << "myidN1 = " << myidN1 << endl;
//       std::cout << "myidN2 = " << myidN2 << endl;
//       std::cout << "myidN3 = " << myidN3 << endl;
      
//     }
    
    if(myidN1 <= myidN2)
    {
      
      All_Proc_Nodes_VP(Node4) = myidN1;
      All_Proc_Nodes_VP(Node4 + Vel_Nnm) = myidN1;
      
      All_Proc_Nodes_Bio(Node4) = myidN1;

    }
    
    if(myidN1 > myidN2)
    {
      
      All_Proc_Nodes_VP(Node4) = myidN2;
      All_Proc_Nodes_VP(Node4 + Vel_Nnm) = myidN2;
      
      All_Proc_Nodes_Bio(Node4) = myidN2;

    }
    
    if(myidN1 <= myidN3)
    {
      
      All_Proc_Nodes_VP(Node6) = myidN1;
      All_Proc_Nodes_VP(Node6 + Vel_Nnm) = myidN1;
      
      All_Proc_Nodes_Bio(Node6) = myidN1;

    }
    
    if(myidN1 > myidN3)
    {
      
      All_Proc_Nodes_VP(Node6) = myidN3;
      All_Proc_Nodes_VP(Node6 + Vel_Nnm) = myidN3;
      
      All_Proc_Nodes_Bio(Node6) = myidN3;

    }
    
    if(myidN2 <= myidN3)
    {
      
      All_Proc_Nodes_VP(Node5) = myidN2;
      All_Proc_Nodes_VP(Node5 + Vel_Nnm) = myidN2;
      
      All_Proc_Nodes_Bio(Node5) = myidN2;

    }
    
    if(myidN2 > myidN3)
    {
      
      All_Proc_Nodes_VP(Node5) = myidN3;
      All_Proc_Nodes_VP(Node5 + Vel_Nnm) = myidN3;
      
      All_Proc_Nodes_Bio(Node5) = myidN3;

    }
  }
  
//   std::cout << "My_Proc_Eles = " << My_Proc_Eles << endl;

  Nodes_Per_Proc_VP = 0;
  Nodes_Per_Proc_Bio = 0;
  
  for(int i = 0; i < Total_Num_Nodes; i++)
  {
    
    if(All_Proc_Nodes_VP(i) == myid)
    {
      
      Nodes_Per_Proc_VP++;
      
    }
  }
  
  for(int i = 0; i < Vel_Nnm; i++)
  {
    
    if(All_Proc_Nodes_Bio(i) == myid)
    {
      
      Nodes_Per_Proc_Bio++;
      
    }
  }
  
  My_Proc_Nodes_VP.Size(Nodes_Per_Proc_VP);
  My_Proc_Nodes_Bio.Size(Nodes_Per_Proc_Bio);
  
  for(int i = 0; i < Total_Num_Nodes; i++)
  {
    
    if(All_Proc_Nodes_VP(i) == myid)
    {
      
      My_Proc_Nodes_VP(new_my_count) = i;
      new_my_count++;
      
    }
    
  }
  
  for(int i = 0; i < Vel_Nnm; i++)
  {
    
    if(All_Proc_Nodes_Bio(i) == myid)
    {
      
      My_Proc_Nodes_Bio(new_my_count_Bio) = i;
      new_my_count_Bio++;
      
    }
    
  }
  
//   if(myid == 3)
//   {
//   
//     std::cout << "All_Proc_Nodes_VP = " << All_Proc_Nodes_VP << endl;
//     std::cout << "Nodes_Per_Proc_VP = " << Nodes_Per_Proc_VP << endl;
//     std::cout << "My_Proc_Nodes_VP = " << My_Proc_Nodes_VP << endl;
//     
//     std::cout << "All_Proc_Nodes_Bio = " << All_Proc_Nodes_Bio << endl;
//     std::cout << "Nodes_Per_Proc_Bio = " << Nodes_Per_Proc_Bio << endl;
//     std::cout << "My_Proc_Nodes_Bio = " << My_Proc_Nodes_Bio << endl;
//     
//     std::cout << "My_Proc_Eles = " << My_Proc_Eles << endl;
//   
//   }
//   
//   int QWERTY;
//   std::cin >> QWERTY;
  
  delete [] Nodes_Per_Proc;
  
}

void ProcNodePartitionClass::ProcNodePartitionVP(int myid,
						 int NumProc,
						 int Vel_Nnm,
						 int Pre_Nnm,
						 int Nem,
						 E_ISDM Vel_Nod, 
						 E_ISDM Pre_Nod,
						 int Vel_Npe,
						 int Pre_Npe,
						 E_ISDM Nod_BC_Hor,
						 E_ISDM Nod_BC_Ver,
						 E_ISDM Pre_Nod_BC)
{
  int *Proc_Nodes_VP = new int[NumProc];

  int *Proc_Eles = new int[NumProc];

  int Total_Num_Nodes, Num_Eles, Mod_Eles;

  Total_Num_Nodes = 2 * Vel_Nnm + Pre_Nnm;

  Mod_Eles = Nem % NumProc;

  Num_Eles = Nem / NumProc; // Since NumProc starts at zero.

  All_Proc_Nodes_VP.Size(Total_Num_Nodes);

  All_Proc_Eles.Size(Nem);

  Shared_Nodes_VP.Size(Total_Num_Nodes);
  
  for(int i = 0; i < Total_Num_Nodes; i++)
  {
    Shared_Nodes_VP(i) = 1.0;
  }
  
  for(int i = 0; i < NumProc; i++)
  {
    if(i < Mod_Eles)
    {
      Proc_Eles[i] = Num_Eles + 1;
    }
    else
    {
      Proc_Eles[i] = Num_Eles;
    }
  }

  if(myid < Mod_Eles)
  {
    Eles_Per_Proc = Num_Eles + 1;
  }
  else
  {
    Eles_Per_Proc = Num_Eles;
  }

  My_Proc_Eles.Size(Eles_Per_Proc);

  int count = 0;

  for(int i = 0; i < NumProc; i++)
  {
    for(int j = 0; j < Proc_Eles[i]; j++)
    {
      if(i == myid) // myid starts at zero.
      {
	All_Proc_Eles(count) = myid;
	My_Proc_Eles(j) = count;
	count++;
      }
      else
      {
	All_Proc_Eles(count) = NumProc + 1;
	count++;
      }
    }
  }

  int Curr_Ele;
  count = 0;

  for(int i = 0; i < Total_Num_Nodes; i++)
  {
    All_Proc_Nodes_VP(i) = NumProc + 1;
  }

  int Curr_Node;

  for(int i = 0; i < Proc_Eles[myid]; i++)
  {
    Curr_Ele = My_Proc_Eles(i);
    
    for(int j = 0; j < Vel_Npe; j++)
    {
      Curr_Node = Vel_Nod(Curr_Ele,j) - 1;
      All_Proc_Nodes_VP(Curr_Node) =  myid;			// Horizontal Velocity Node
      All_Proc_Nodes_VP(Curr_Node + Vel_Nnm) =  myid;		// Vertical Velocity Node
      
      count = count + 2;
      
    }
    
    for(int j = 0; j < Pre_Npe; j++)
    {
      Curr_Node = Pre_Nod(Curr_Ele,j) - 1;
      All_Proc_Nodes_VP(Curr_Node + 2 * Vel_Nnm) =  myid;	// Pressure Node
      
      count++;
      
    }
  }
  
  for(int e = 0; e < Nem; e++)
  {
    if(All_Proc_Eles(e) != myid)
    {
      for(int i = 0; i < Vel_Npe; i++)
      {
	if(All_Proc_Nodes_VP(Vel_Nod(e,i) - 1) == myid)
	{
	  Shared_Nodes_VP(Vel_Nod(e,i) - 1) = 2.0;		// Horizontal Velocity
	  Shared_Nodes_VP(Vel_Nod(e,i) - 1 + Vel_Nnm) = 2.0;	// Vertical Velocity
	}
      }
    }
  }

//   std::cout << "myid = " << myid << endl;
//   std::cout << " Shared_Nodes_VP = " << Shared_Nodes_VP << endl;
//   
//   int QWERT;
//   std::cin >> QWERT;
  
  count = 0;

  for(int i = 0; i < Total_Num_Nodes; i++)
  {
    if(All_Proc_Nodes_VP(i) == myid)
    {
      count++;
    }
  }
  
  My_Proc_Nodes_VP.Size(count);
  
  Nodes_Per_Proc_VP = count;
  
  count = 0;
  
  for(int i = 0; i < Total_Num_Nodes; i++)
  {
    if(All_Proc_Nodes_VP(i) == myid)
    {
      
      My_Proc_Nodes_VP(count) = i;
      count++;
      
    }
  }

  if(myid == 1)
  {
    std::cout << "myid = " << myid << endl;
    std::cout << "All_Proc_Nodes = " << All_Proc_Nodes_VP << endl;
    std::cout << "Nodes_Per_Proc = " << Nodes_Per_Proc_VP << endl;
    std::cout << "My_Proc_Nodes = " << My_Proc_Nodes_VP << endl;
    
    std::cout << "All_Proc_Eles = " << All_Proc_Eles << endl;
    std::cout << "Eles_Per_Proc = " << Eles_Per_Proc << endl;
    std::cout << "My_Proc_Eles = " << My_Proc_Eles << endl;
    
    int HKJL;
    std::cin >> HKJL;
  }
  
  delete [] Proc_Nodes_VP;
  delete [] Proc_Eles;
  
}

void ProcNodePartitionClass::ProcNodePartitionBio(int myid, 
						 int NumProc, 
						 int Vel_Nnm, 
						 int Pre_Nnm, 
						 int Nem,
						 E_ISDM Vel_Nod, 
						 E_ISDM Pre_Nod,
						 int Vel_Npe,
						 int Pre_Npe,
						 E_ISDM Bio_Nod_BC)
{
  int *Proc_Nodes_Bio = new int[NumProc];
  
  int *Proc_Eles = new int[NumProc];
 
  int Total_Num_Nodes, Num_Eles, Mod_Eles;

  Total_Num_Nodes = Vel_Nnm;

  Mod_Eles = Nem % NumProc;
  
  Num_Eles = Nem / NumProc; // Since NumProc starts at zero.

  All_Proc_Nodes_Bio.Size(Total_Num_Nodes);
  
  All_Proc_Eles.Size(Nem);
  
  Shared_Nodes_Bio.Size(Total_Num_Nodes);
  
  for(int i = 0; i < Total_Num_Nodes; i++)
  {
    Shared_Nodes_Bio(i) = 1.0;
  }

  for(int i = 0; i < NumProc; i++)
  {
    if(i < Mod_Eles)
    {
      Proc_Eles[i] = Num_Eles + 1;
    }
    else
    {
      Proc_Eles[i] = Num_Eles;
    }
  }

  if(myid < Mod_Eles)
  {
    Eles_Per_Proc = Num_Eles + 1;
  }
  else
  {
    Eles_Per_Proc = Num_Eles;
  }
  
  My_Proc_Eles.Size(Eles_Per_Proc);
  
  int count = 0;
  
  for(int i = 0; i < NumProc; i++)
  {
    for(int j = 0; j < Proc_Eles[i]; j++)
    {
      if(i == myid) // myid starts at zero.
      {
	All_Proc_Eles(count) = myid;
	My_Proc_Eles(j) = count;
	count++;
      }
      else
      {
	All_Proc_Eles(count) = NumProc + 1;
	count++;
      }
    }
  }
  
//   if(myid == 0)
//   {
//     std::cout << "myid = " << myid << endl;
//     std::cout << "All_Proc_Eles = " << All_Proc_Eles << endl;
//     std::cout << "Eles_Per_Proc = " << Eles_Per_Proc << endl;
//     std::cout << "My_Proc_Eles = " << My_Proc_Eles << endl;
//   }
//   
//   std::cout << "Here 1" << endl;
  
  int Curr_Ele;
  count = 0;
  
  for(int i = 0; i < Total_Num_Nodes; i++)
  {
    All_Proc_Nodes_Bio(i) = NumProc + 1;
  }
  
  int Curr_Node;
  
  for(int i = 0; i < Proc_Eles[myid]; i++)
  {
    Curr_Ele = My_Proc_Eles(i);
    
    for(int j = 0; j < Vel_Npe; j++)
    {
      Curr_Node = Vel_Nod(Curr_Ele,j) - 1;
      All_Proc_Nodes_Bio(Curr_Node) =  myid;
      
      count++;
      
    }
  }

  for(int e = 0; e < Nem; e++)
  {
//     std::cout << "e = " << e << endl;
    
    if(All_Proc_Eles(e) != myid)
    {
      for(int i = 0; i < Vel_Npe; i++)
      {
// 	std::cout << "i = " << i << endl;
	
	if(All_Proc_Nodes_Bio(Vel_Nod(e,i) - 1) == myid)
	{
	  Shared_Nodes_Bio(Vel_Nod(e,i) - 1) = 2.0;
	}
      }
    }
  }
  
//   std::cout << "myid = " << myid << endl;
//   std::cout << " Shared_Nodes_Bio = " << Shared_Nodes_Bio << endl;
//   
//   int QWERT;
//   std::cin >> QWERT;
  
  count = 0;

  for(int i = 0; i < Total_Num_Nodes; i++)
  {
    if(All_Proc_Nodes_Bio(i) == myid)
    {
      count++;
    }
  }
  
  My_Proc_Nodes_Bio.Size(count);
  
  Nodes_Per_Proc_Bio = count;
  
  count = 0;
  
  for(int i = 0; i < Total_Num_Nodes; i++)
  {
    if(All_Proc_Nodes_Bio(i) == myid)
    {
      
      My_Proc_Nodes_Bio(count) = i;
      count++;
      
    }
  }
  
//   std::cout << "Here 3" << endl;
//   
//   if(myid == 0)
//   {
//     std::cout << "myid = " << myid << endl;
//     std::cout << "All_Proc_Nodes = " << All_Proc_Nodes_Bio << endl;
//     std::cout << "Nodes_Per_Proc = " << Nodes_Per_Proc_Bio << endl;
//     std::cout << "My_Proc_Nodes = " << My_Proc_Nodes_Bio << endl;
//     
//     int HKJL;
//     std::cin >> HKJL;
//   }
  
  delete [] Proc_Nodes_Bio;
  delete [] Proc_Eles;
}