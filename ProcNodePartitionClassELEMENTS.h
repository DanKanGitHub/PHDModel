#ifndef PROCNODEPARTITIONCLASS_H_
#define PROCNODEPARTITIONCLASS_H_

#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseMatrix.h"

typedef Epetra_IntSerialDenseVector E_ISDV;
typedef Epetra_SerialDenseVector E_SDV;
typedef Epetra_IntSerialDenseMatrix E_ISDM;

class ProcNodePartitionClass
{
  
private:
  
public:
  ProcNodePartitionClass();
  ~ProcNodePartitionClass();

  void ProcNodePartitionVP(int myid, 
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
			   E_ISDM Pre_Nod_BC);
  
  E_ISDV All_Proc_Nodes_VP;
  int Nodes_Per_Proc_VP;
  E_ISDV My_Proc_Nodes_VP;
  
  E_ISDV All_Proc_Eles;
  int Eles_Per_Proc;
  E_ISDV My_Proc_Eles;
  
  void ProcNodePartitionBio(int myid, 
			   int NumProc, 
			   int Vel_Nnm, 
			   int Pre_Nnm, 
			   int Nem, 
			   E_ISDM Vel_Nod, 
			   E_ISDM Pre_Nod, 
			   int Vel_Npe,
			   int Pre_Npe,
			   E_ISDM Bio_Nod_BC);
  
  E_ISDV All_Proc_Nodes_Bio;
  int Nodes_Per_Proc_Bio;
  E_ISDV My_Proc_Nodes_Bio;
  
  E_SDV Shared_Nodes_VP, Shared_Nodes_Bio;
};

#endif