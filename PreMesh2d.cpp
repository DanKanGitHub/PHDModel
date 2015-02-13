// mesh generation function of a rectangular domain for 2D finite element
// method
// Input: XL, XR (x-coordinates of left and right boundary)
//        YB, YT (y-coordinates of bottom and top boundary)
//        NX, NY (number of intervals in x and y direction)
//        FLAG = 1: linear triangle,
//		 2: quadratic triangle.
//	  Nem: # of elements in the mesh
// Output: Nnm: # of nodes in the mesh
//         Npe: # of nodes per element
//         Glxy: coordinates of global nodes, array of size Nnmx2
//         	Glxy(I,0) = x-coord of I-th global node, Glxy(I,1) = y-coord
//         Nod: connectivity matrix, array of size NemxNpe
//         Nb_Ele_No: number of elements with Natural BC
//         Nb_Ele: element list with natural BC, array of size Nb_Ele_Nox2
//         	   Nb_Ele(I,1): element index of I-th Natural boundary element
//         	   Nb_Ele(I,2): number of natural boundary condition edges with NBC in I-th boundary element
//         Nb_Edg: local edge number of the boundary element, array of size
//         	   Nb_Ele_Nox2, NB_EDG(I,j) = local edge number of j-th
//                 boundary edge in I-the boundary element, 1 <= j <= NBELE(I,2)
//         Eb_Node_No: number of nodes with Essential BCs
//         Eb_Node: global node number of EBC node, array of size EB_NODE_NOx1

// Element and node counting starts at 1.
// Example:
// |'''/|'''/|
// | 6/ | 8/ |
// | / 5| / 7|
// |/___|/___|
// |   /|   /|
// | 2/ | 4/ |
// | / 1| / 3|
// |/___|/___|

// Quadrtaic triangle node structure
//              3
//             /|
//            / |
//           /  |
// Eta < 0 6/   |5 Xi < 0
//         /    |
//        /     |
//       /______|
//       1  4   2
//     Xi + Eta > 1

//           Xi < 0
//          3  5   2
//          |'''''/
//          |    / 
//          |   /  
// Eta < 0 6|  /4 Xi + Eta > 1
//          | /    
//          |/  
//          /
//          1

// Barycentric Triangle
//
//     2
//     |\
//     | \
// Eta |  \
//     |   \
//     |    \
//    3|_____\1
//        Xi


// Boundary Condition Flag:
// 0 means None
// 1 means Essential
// 2 means Natural

#include "PreMesh2d.h"

void PreMesh2d(double dx, 
		double dy, 
		double XL,
		double YB,
		int NX, 
		int NY,
		int Nem, 
		int & Nnm, 
		int & Npe, 
		E_SDM & Glxy, 
		E_ISDM & Nod,
		E_ISDM & Nod_BC) {
  // Variable declarations
  int Eln1, Eln2;
  
  // % size of various structures
  Nnm = (NY + 1) * (NX + 1);
  Npe = 3;
  
  Glxy.Shape(Nnm,2);
  for(int i = 0; i <= Nnm - 1; i++)
  {
    Glxy(i,0) = 0.0;
    Glxy(i,1) = 0.0;
  }
  
  Nod.Shape(Nem,Npe); // Added a dimension for the EBC, NBC flag.
  Nod_BC.Shape(Nem,Npe); //
  // Initialize
//     for(int i = 0; i <= Nem - 1; i++)
//     {
//       for(int j = 0; j <= Npe - 1; j++)
//       {
// 	Nod(i,j,0) = 0;
//       }
//     }
  
  // global coordinate of nodes
  for (int i = 0; i <= Nnm-1; i++)
  {
    Glxy(i,0) = (i % (NX + 1)) * dx + XL;
    Glxy(i,1) = floor(i / (NX + 1)) * dy + YB;
  }

  // connectivity matrix
  for (int j = 0; j <= NY-1; j++)	// indexes the horizontal elements
  {
    for (int i = 0; i <= NX-1; i++)	// indexes the vertical elements
    {
      Eln1 = j*2*NX + 2*(i + 1) - 1;	// lower element number
      Eln2 = Eln1 + 1;		// upper element number

      Nod(Eln1-1,0) = j*(NX+1) + i + 1;
      Nod(Eln1-1,1) = Nod(Eln1-1,0) + 1;
      Nod(Eln1-1,2) = Nod(Eln1-1,1) + NX + 1;

      Nod(Eln2-1,0) = Nod(Eln1-1,0);
      Nod(Eln2-1,1) = Nod(Eln1-1,2);
      Nod(Eln2-1,2) = Nod(Eln2-1,1) - 1;
    }
  }
  
  // There is only one pressure BC, it is essential, and at the first global node
//   Nod_BC(0,0) = 1; // Lower element local node one.
//   Nod_BC(1,0) = 1; // Upper element local node one.

}