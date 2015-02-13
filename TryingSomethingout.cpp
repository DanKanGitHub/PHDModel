// Processor node Partition

#include "ProcNodePartitionClass.h"

ProcNodePartitionClass::ProcNodePartitionClass()
{

}

ProcNodePartitionClass::~ProcNodePartitionClass()
{

}

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

  int Total_Num_Nodes, Num_Nodes, Mod_Nodes, count;

  Total_Num_Nodes = 2 * Vel_Nnm + Pre_Nnm;
  
  // Partition the pressure nodes first.
  Mod_Nodes = Pre_Nnm % NumProc;
  
  Num_Nodes = Pre_Nnm / NumProc; // Since NumProc starts at zero.
  
  int *Nodes_Per_Proc = new int[NumProc];
  int *Pre_Node_Part = new int[Pre_Nnm];
  
  All_Proc_Nodes_VP.Size(Total_Num_Nodes);
  All_Proc_Nodes_Bio.Size(Vel_Nnm);
  
  // Pressure Num of nodes partitiioned
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
  
  // Assign a check value  Might be unnecessary
//   for(int i = 0; i < Total_Num_Nodes; i++)
//   {
//     All_Proc_Nodes_VP(i) = NumProc;
//   }
//   
//   for(int i = 0; i < Vel_Nnm; i++)
//   {
//     All_Proc_Nodes_Bio(i) = NumProc;
//   }
  
  count = 0;
  
  // Assign pressure nodes myids
  for(int i = 0; i < NumProc; i++)
  {
    for(int j = 0; j < Nodes_Per_Proc[i]; j++)
    {
    
      Pre_Node_Part(count) = i;
      All_Proc_Nodes_VP(count + 2 * Vel_Nnm) = i;
      count++;
      
    }
  }
  
  // Find associated velocity nodes for previously assigned pressure nodes
  
  for(int i = 0; i < NumProc; i++)
  {
    for(int e = 0; e < Nem; e++)
    {
      for(int n = 0; n < Pre_Npe; n++)
      {
	if(Pre_Node_Part(Pre_Nod(e,n) - 1) == i)
	{
	  
	  All_Proc_Nodes_VP(Vel_Nod(e,n) - 1) = i;
	  All_Proc_Nodes_VP(Vel_Nod(e,n) - 1 + Vel_Nnm) = i;
	  
	  All_Proc_Nodes_Bio(Vel_Nod(e,n) - 1) = i;
	  
	}
      }
    }
  }

  int Node1, Node2, Node3, Node4, Node5, Node6, myidN1, myidN2, myidN3;
  
  // Assign non-vertex nodes
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
  
  Nodes_Per_Proc_VP = 0;
  Nodes_Per_Proc_Bio = 0;
  
  // Count how many nodes "I" own
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
  
  // Build local global node arrays
  My_Proc_Nodes_VP.Size(Nodes_Per_Proc_VP);
  My_Proc_Nodes_Bio.Size(Nodes_Per_Proc_Bio);
  
  int new_my_count, new_my_count_Bio;
  
  new_my_count = 0;
  new_my_count_Bio = 0;
  
  //Assign local nodes to arrays
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
  
  delete [] Nodes_Per_Proc;
  delete [] Pre_Node_Part;
}