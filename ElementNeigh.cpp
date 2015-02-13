// ElementNeigh determines all the elements which share a node
// The numbering is ddetermined by the local node number, i.e., side 1 is opposite node 1.
// This only works for my mesh
// NX is the number of intervals in the horizontal direction.

// Quadrtaic triangle node structure
// Odd Element
//        3
//       /|
//      / |
//     /  |
//   6/   |5
//   /    |
//  /     |
// /______|
// 1  4   2


// Even Element
//  3  5   2
//  |'''''/
//  |    / 
//  |   /  
// 6|  /4
//  | /    
//  |/  
//  1

// counting of elements starts at 1.
// The side opposite the nodes 1,2,3 are entered into Ele_Neigh in the 1st, 2nd and 3rd columns respectively.

#include "ElementNeigh.h"
#include "Epetra_IntSerialDenseMatrix.h"

typedef Epetra_IntSerialDenseMatrix E_ISDM;

void ElementNeigh(int Nem,
		  int NX,
		  E_ISDM & Ele_Neigh)
{
  NX = 2*NX; // *2 because we have 2 triangles in each interval
  
  for (int i = 1; i <= Nem; i++)
  {
    if(i % 2 > 0) // i.e., odd
    {
      // if element i is at the end of a row then (i+1) % NX = 0 otherwise it is postive.
      if((i+1) % NX > 0) 
      {
	Ele_Neigh(i-1,0) = i + 3;
      }
      else // if this neighbor does not exsist
      {
	Ele_Neigh(i-1,0) = -1; // This acts as a check for later use
      }
      
      // Every odd element has this neighbor
      Ele_Neigh(i-1,1) = i + 1;
      
      // If the element is not in the first row.
      if(i-NX+1 > 0) 
      {
	Ele_Neigh(i-1,2) = i - NX + 1;
      }
      else
      {
	Ele_Neigh(i-1,2) = -1; // 1st row doesn't have this neighbor
      }
    }
    else // element is even
    {
      // if element i is at the start of a row then (i-2) % NX = 0 otherwise it is postive.
      if((i-2) % NX > 0) 
      {
	Ele_Neigh(i-1,1) = i - 3;
      }
      else
      {
	Ele_Neigh(i-1,1) = -1;
      }
      
      // Every even element has this neighbor
      Ele_Neigh(i-1,2) = i - 1;
      
      // If element is not in the last row.
      if(i+NX-1 < Nem) 
      {
	Ele_Neigh(i-1,0) = i + NX - 1;
      }
      else
      {
	Ele_Neigh(i-1,0) = -1; // Last row doesn't have this neighbor
      } // if
    } // if odd
  } // for
} // function