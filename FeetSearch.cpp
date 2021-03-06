// Determines which element the point x is in.
// Cur_Ele:	Current element number
// New_Ele:	New element number
// Ele_Neigh:	All elements which are the Neighbor of CEle, size 4XNEM

// Element counting starts at 1 not zero

#include "FeetSearch.h"
#include "Epetra_IntSerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

typedef Epetra_SerialDenseVector E_SDV;
typedef Epetra_IntSerialDenseMatrix E_ISDM;

int FeetSearch(E_SDV x, 
	       int Cur_Ele, 
	       E_ISDM Ele_Neigh)
{
  int New_Ele;
  
  New_Ele = Cur_Ele;
  
  if(Cur_Ele % 2 > 0) // i.e., odd
  {
    if(x(0) < 0)
    {
      New_Ele = Ele_Neigh(Cur_Ele - 1, 0); // -1 for C++ indexing
    }
    else if(x(0) + x(1) > 1)
    {
      New_Ele = Ele_Neigh(Cur_Ele - 1, 2);
    }
    else if(x(1) < 0)
    {
      New_Ele = Ele_Neigh(Cur_Ele - 1, 1);
    }
  }
  else
  {
    if(x(0) < 0)
    {
      New_Ele = Ele_Neigh(Cur_Ele - 1, 0); // -1 for C++ indexing
    }
    else if(x(0) + x(1) > 1)
    {
      New_Ele = Ele_Neigh(Cur_Ele - 1, 2);
    }
    else if(x(1) < 0)
    {
      New_Ele = Ele_Neigh(Cur_Ele - 1, 1);
    }
  }
 
  return New_Ele;
}