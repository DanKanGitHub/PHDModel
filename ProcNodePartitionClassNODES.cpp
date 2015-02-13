// Processor node Partition

#include "ProcNodePartitionClass.h"

ProcNodePartitionClass::ProcNodePartitionClass()
{
  
  std::cout << "Constructor" << endl;
  
}

ProcNodePartitionClass::~ProcNodePartitionClass()
{
  
  std::cout << "Destructor" << endl;
  
}

void ProcNodePartitionClass::ProcNodePartition(int myid, int NumProc, int Vel_Nnm, int Pre_Nnm)
{
 
  int Total_Num_Nodes, Num_Nodes, Mod_Nodes;

  Total_Num_Nodes = 2 * Vel_Nnm + Pre_Nnm;
  
  Mod_Nodes = Total_Num_Nodes % (NumProc + 1);
  
  Num_Nodes = Total_Num_Nodes / (NumProc + 1); // Since NumProc starts at zero.
  
  int *Nodes_Per_Proc = new int[NumProc + 1];
  
  All_Proc_Nodes.Size(Total_Num_Nodes);
  
  for(int i = 0; i <= NumProc; i++)
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
  
  int count = 0;
  
  for(int i = 0; i < NumProc; i++)
  {
    for(int j = 0; j < Nodes_Per_Proc[i]; j++)
    {
      if(i == myid)
      {
	All_Proc_Nodes(count) = myid;
	count++;
      }
      else
      {
	All_Proc_Nodes(count) = 0;
	count++;
      }
    }
  }
  
  delete [] Nodes_Per_Proc;
  
}