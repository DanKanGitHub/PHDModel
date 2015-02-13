// Processor node Partition

#include "ProcNodePartitionClass.h"

ProcNodePartitionClass::ProcNodePartitionClass()
{

}

ProcNodePartitionClass::~ProcNodePartitionClass()
{

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

//   if(myid == 0)
//   {
//     std::cout << "myid = " << myid << endl;
//     std::cout << "All_Proc_Nodes = " << All_Proc_Nodes_VP << endl;
//     std::cout << "Nodes_Per_Proc = " << Nodes_Per_Proc_VP << endl;
//     std::cout << "My_Proc_Nodes = " << My_Proc_Nodes_VP << endl;
//     
//     std::cout << "All_Proc_Eles = " << All_Proc_Eles << endl;
//     std::cout << "Eles_Per_Proc = " << Eles_Per_Proc << endl;
//     std::cout << "My_Proc_Eles = " << My_Proc_Eles << endl;
//     
//     int HKJL;
//     std::cin >> HKJL;
//   }
  
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