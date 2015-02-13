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

#include "VelBioMesh2d.h"

void VelBioMesh2d(double dx, 
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
		  E_ISDM & Nod_BC_Hor,
		  E_ISDM & Nod_BC_Ver,
		  E_ISDM & Bio_Nod_BC) {
    // Variable declarations
  int Eln1, Eln2;
  
  NX = 2 * NX;
  NY = 2 * NY;
  
  dx = dx / 2;
  dy = dy / 2;
  
  // size of various structures
  Nnm = (NY + 1) * (NX + 1);
  Npe = 6;
  
  Glxy.Shape(Nnm,2);
  for(int i = 0; i <= Nnm - 1; i++)
  {
    Glxy(i,0) = 0.0;
    Glxy(i,1) = 0.0;
  }
  
  Nod.Shape(Nem,Npe);
  Nod_BC_Hor.Shape(Nem,Npe);
  Nod_BC_Ver.Shape(Nem,Npe);
  Bio_Nod_BC.Shape(Nem,Npe);
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
  for (int j = 0; j <= NY/2-1; j++) 	// indexes the horizontal elements
  {
    for (int i = 0; i <= NX/2-1; i++)// indexes the vertical elements
    {
      Eln1 = j*NX + 2*(i + 1) - 1;	// lower element number,
      Eln2 = Eln1 + 1;		// upper element number

      Nod(Eln1-1,0) = 2 * j * (NX + 1) + 2 * i + 1; //Global node number
      Nod(Eln1-1,1) = Nod(Eln1-1,0) + 2;
      Nod(Eln1-1,2) = Nod(Eln1-1,1) + 2*(NX + 1);
      Nod(Eln1-1,3) = Nod(Eln1-1,0) + 1;
      Nod(Eln1-1,4) = Nod(Eln1-1,1) + NX + 1;
      Nod(Eln1-1,5) = Nod(Eln1-1,4) - 1;

      Nod(Eln2-1,0) = Nod(Eln1-1,0);
      Nod(Eln2-1,1) = Nod(Eln1-1,2);
      Nod(Eln2-1,2) = Nod(Eln2-1,1) - 2;
      Nod(Eln2-1,3) = Nod(Eln1-1,5);
      Nod(Eln2-1,4) = Nod(Eln2-1,1) - 1;
      Nod(Eln2-1,5) = Nod(Eln2-1,3) - 1;
      
      // Specific to the mesh I generated
      // Node one from the lower element On the bottom row
      if( ( Nod(Eln1-1,0) >= 1) & (Nod(Eln1-1,0) <= NX) )
      {
	//Lower Element
	Nod_BC_Hor(Eln1-1,0) = 1; // Flag = 1 means Essential BC
	Nod_BC_Hor(Eln1-1,1) = 1;
	Nod_BC_Hor(Eln1-1,3) = 1;
	
	Bio_Nod_BC(Eln1-1,0) = 2; // Flag = 1 means Essential BC
	Bio_Nod_BC(Eln1-1,1) = 2;
	Bio_Nod_BC(Eln1-1,3) = 2;
	
	//Upper Element
	Nod_BC_Hor(Eln2-1,0) = 1;
	
	Bio_Nod_BC(Eln2-1,0) = 2;
	
	//Lower Element
	Nod_BC_Ver(Eln1-1,0) = 1; // Flag = 1 means Essential BC
	Nod_BC_Ver(Eln1-1,1) = 1;
	Nod_BC_Ver(Eln1-1,3) = 1;
	
	//Upper Element
	Nod_BC_Ver(Eln2-1,0) = 1;
      }

      // Nod 3 from the lower element On the Top row
      if( ( Nod(Eln1-1,2) >= NX * (NY + 1) ) & ( Nod(Eln1-1,2) <= (NX + 1) * (NY + 1) ) )
      {
	//Lower Element
	Nod_BC_Hor(Eln1-1,2) = 1; // Flag = 1 means Essential BC
	
	Bio_Nod_BC(Eln1-1,2) = 2; // Flag = 1 means Essential BC
	
	//Upper Element
	Nod_BC_Hor(Eln2-1,1) = 1;
	Nod_BC_Hor(Eln2-1,2) = 1;
	Nod_BC_Hor(Eln2-1,4) = 1;
	
	Bio_Nod_BC(Eln2-1,1) = 2;
	Bio_Nod_BC(Eln2-1,2) = 2;
	Bio_Nod_BC(Eln2-1,4) = 2;
	
	//Lower Element
	Nod_BC_Ver(Eln1-1,2) = 1; // Flag = 1 means Essential BC
	
	//Upper Element
	Nod_BC_Ver(Eln2-1,1) = 1;
	Nod_BC_Ver(Eln2-1,2) = 1;
	Nod_BC_Ver(Eln2-1,4) = 1;
	
      }
      
      // Node one from the lower element On the front
      // By doing the front last it ensures that the Bio_BC on the bottom and top front corners
      // are labeled as essential BCs
      if( (Nod(Eln1-1,0)) % (NX + 1) == 1 )
      {
	//Lower Element
	Nod_BC_Hor(Eln1-1,0) = 1; // Flag = 1 means Essential BC
	
	Bio_Nod_BC(Eln1-1,0) = 1; // Flag = 1 means Essential BC
	
	//Upper Element
	Nod_BC_Hor(Eln2-1,0) = 1;
	Nod_BC_Hor(Eln2-1,2) = 1;
	Nod_BC_Hor(Eln2-1,5) = 1;
	
	Bio_Nod_BC(Eln2-1,0) = 1;
	Bio_Nod_BC(Eln2-1,2) = 1;
	Bio_Nod_BC(Eln2-1,5) = 1;
	
	//Lower Element
	Nod_BC_Ver(Eln1-1,0) = 1; // Flag = 1 means Essential BC
	
	//Upper Element
	Nod_BC_Ver(Eln2-1,0) = 1;
	Nod_BC_Ver(Eln2-1,2) = 1;
	Nod_BC_Ver(Eln2-1,5) = 1;
      }
      
      // Nod 3 from the lower element On the back above and below the bottom and top rows respectively
      if( ( Nod(Eln1-1,2) > (NX + 1) ) & ( Nod(Eln1-1,2) <= (NX + 1) * (NY + 1) ) & ( Nod(Eln1-1,2) % (NX + 1) == 0) )
      {
	if(Eln1 == (NX - 1))
	{
	  //Lower Element
	  Nod_BC_Hor(Eln1-1,2) = 2; // Flag = 2 means Natural BC
	  Nod_BC_Hor(Eln1-1,4) = 2;
	  
	  Bio_Nod_BC(Eln1-1,2) = 2; // Flag = 2 means Natural BC
	  Bio_Nod_BC(Eln1-1,4) = 2;
	  
	  //Upper Element
	  Nod_BC_Hor(Eln2-1,1) = 2;
	  
	  Bio_Nod_BC(Eln2-1,1) = 2;
	  
	  //Lower Element
	  Nod_BC_Ver(Eln1-1,2) = 1; // Flag = 2 means Natural BC
	  Nod_BC_Ver(Eln1-1,4) = 1;
	  
	  //Upper Element
	  Nod_BC_Ver(Eln2-1,1) = 1;
	}
	else if(Eln1 == (NX * (NY / 2) - 1))
	{
	  //Lower Element
	  Nod_BC_Hor(Eln1-1,1) = 2;
	  Nod_BC_Hor(Eln1-1,4) = 2;
	  
	  Bio_Nod_BC(Eln1-1,1) = 2;
	  Bio_Nod_BC(Eln1-1,4) = 2;
	  
	  //Lower Element
	  Nod_BC_Ver(Eln1-1,1) = 1;
	  Nod_BC_Ver(Eln1-1,4) = 1;
	}
	else
	{
	  //Lower Element
	  Nod_BC_Hor(Eln1-1,1) = 2;
	  Nod_BC_Hor(Eln1-1,2) = 2; // Flag = 2 means Natural BC
	  Nod_BC_Hor(Eln1-1,4) = 2;
	  
	  Bio_Nod_BC(Eln1-1,1) = 2;
	  Bio_Nod_BC(Eln1-1,2) = 2; // Flag = 2 means Natural BC
	  Bio_Nod_BC(Eln1-1,4) = 2;
	  
	  //Upper Element
	  Nod_BC_Hor(Eln2-1,1) = 2;
	  
	  Bio_Nod_BC(Eln2-1,1) = 2;
	  
	  //Lower Element
	  Nod_BC_Ver(Eln1-1,1) = 1;
	  Nod_BC_Ver(Eln1-1,2) = 1; // Flag = 2 means Natural BC
	  Nod_BC_Ver(Eln1-1,4) = 1;
	  
	  //Upper Element
	  Nod_BC_Ver(Eln2-1,1) = 1;
	}
      }
    }
  }
}