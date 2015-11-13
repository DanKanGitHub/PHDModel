// Processor node Partition

#include "ProcNodePartitionClass.h"

ProcNodePartitionClass::ProcNodePartitionClass()
{

}

ProcNodePartitionClass::~ProcNodePartitionClass()
{

}

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

  int Total_Num_Nodes, Num_Nodes, Mod_Nodes, EleNum;
  int count, NpeNum, Global_Node_VP, Global_Node_Pre;
  int Node1, Node2, Node3, Node4, Node5, Node6, myidN1, myidN2, myidN3;
  int new_my_count, new_my_count_Bio;

  Total_Num_Nodes = 2 * Vel_Nnm + Pre_Nnm;
  
  // Partition the pressure nodes first because the pressure nodes are the vertices
  // of each element and the processor assignment on the vertices govern the rest of the nodes.
  Mod_Nodes = Pre_Nnm % NumProc; // remainder
  
  Num_Nodes = Pre_Nnm / NumProc; // Integer division
  
  int *Nodes_Per_Proc = new int[NumProc];
  
  All_Proc_Nodes_VP.Size(Total_Num_Nodes);
  All_Proc_Nodes_Bio.Size(Vel_Nnm);
  
  for(int i = 0; i < NumProc; i++)
  {
    if(i < Mod_Nodes)
    {
      Nodes_Per_Proc[i] = Num_Nodes + 1; // for the number of processors that make up the remainder
    }
    else
    {
      Nodes_Per_Proc[i] = Num_Nodes;
    }
  }

  // Initialize to a value that myid will never take
  for(int i = 0; i < Total_Num_Nodes; i++)
  {
    All_Proc_Nodes_VP(i) = NumProc;
  }
  
  for(int i = 0; i < Vel_Nnm; i++)
  {
    All_Proc_Nodes_Bio(i) = NumProc;
  }

  // Initialize counters
  EleNum = 0;
  NpeNum = 0;
  
  // Loop over processors for node assignment to processors
  for(int i = 0; i < NumProc; i++)
  {
    
    count = 0;

    // Assign nodes to a processor until count is less than the number of nodes
    // we determined the processor is to have. We use strictly less than because 
    // count starts at zero and Nodes_Per_Proc starts counting at 1
    while(count < Nodes_Per_Proc[i])
    {
      
      // Once we have looped through all the  nodes in an element then index the 
      // element the number
      if(NpeNum == Pre_Npe)
      {
	
	NpeNum = 0;
	EleNum++;
	
      }

      Global_Node_VP = Vel_Nod(EleNum,NpeNum) - 1;
      Global_Node_Pre = Pre_Nod(EleNum,NpeNum) - 1;

      // If a node hasn't been assigned yet.
      if(All_Proc_Nodes_VP(Global_Node_VP) == NumProc)
      {
	All_Proc_Nodes_VP(Global_Node_Pre + 2 * Vel_Nnm) = i;	// Pressure comes after all the vel nodes
	All_Proc_Nodes_VP(Global_Node_VP) = i;			// Hor Vel
	All_Proc_Nodes_VP(Global_Node_VP + Vel_Nnm) = i;	// Vert Vel comes after all the Hor Vel nodes
	
	All_Proc_Nodes_Bio(Global_Node_VP) = i;

	count++; // Only count when a node is assigned
      }
      
      NpeNum++;
      
    } // while
  } // for

  My_Proc_Eles.Size(Nem);
  
  // Initialize to a value that myid will never equal
  for(int i = 0; i < Nem; i++)
  {
    
    My_Proc_Eles(i) = -1;
    
  }

  // Now do the nodes of each element that are not the vertices
  for(int i = 0; i < Nem; i++)
  {
    
    // Global node numbers for the nodes of the element
    Node1 = Vel_Nod(i,0) - 1;
    Node2 = Vel_Nod(i,1) - 1;
    Node3 = Vel_Nod(i,2) - 1;
    Node4 = Vel_Nod(i,3) - 1;
    Node5 = Vel_Nod(i,4) - 1;
    Node6 = Vel_Nod(i,5) - 1;
    
    // Processor ids for the vertices
    myidN1 = All_Proc_Nodes_VP(Node1);
    myidN2 = All_Proc_Nodes_VP(Node2);
    myidN3 = All_Proc_Nodes_VP(Node3);
    
    // Assign a processor an element
    if(myid == myidN1 || myid == myidN2 || myid == myidN3)
    {
      
      My_Proc_Eles(i) = myid;
      
    }

    // Assign non-vertices nodes to a processor with the lowest
    // id number.
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
  
  // Count how many nodes each processor owns
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
  
  new_my_count = 0;
  new_my_count_Bio = 0;
  
  // Assign global node numbers for each node owned by a processor
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
  
}