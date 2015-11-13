#include "StressReassemble.h"

void StressReassemble(int Str_Nnm,
		      int myid,
		      E_SDM Str_New, 
		      E_ISDV All_Proc_Nodes, 
		      E_SDM & Str)
{
  for(int i = 0; i < Str_Nnm; i++)
  {
    if(All_Proc_Nodes(i) == myid)
    {
      Str(i,0) = Str_New(i,0);
      Str(i,1) = Str_New(i,1);
      Str(i,2) = Str_New(i,2);
    }
  }
}