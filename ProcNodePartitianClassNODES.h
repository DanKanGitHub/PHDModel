#ifndef PROCNODEPARTITIONCLASS_H_
#define PROCNODEPARTITIONCLASS_H_

#include "Epetra_IntSerialDenseVector.h"

class ProcNodePartitionClass
{
  
private:
  
public:
  ProcNodePartitionClass();
  ~ProcNodePartitionClass();
  void ProcNodePartition(int myid, int NumProc, int Vel_Nnm, int Pre_Nnm);
  Epetra_IntSerialDenseVector All_Proc_Nodes;
  
}

#endif
