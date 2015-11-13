// mesh generation function of a rectangular domain for 2D finite element
// method
// Input: XL, XR (x-coordinates of left and right boundary)
//        YB, YT (y-coordinates of bottom and top boundary)
//        NX, NY (number of intervals in x and y direction)
//        FLAG = 1: linear triangle,
//		 2: quadratic triangle.
//	  Nem: # of elements in the mesh
// Output: TotalNumNod: # of nodes in the mesh
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
//               3
//              /|
//             / |
//            /  |
// Eta < 0, 6/   |5, Xi < 0
//          /    |
//         /     |
//        /______|
//        1  4   2
//     Xi + Eta > 1

//           Xi < 0
//           3  5   2
//           |'''''/
//           |    / 
//           |   /  
// Eta < 0, 6|  /4, Xi + Eta > 1
//           | /    
//           |/  
//           /
//           1

// Barycentric Triangle
//
//      2
//      |\
//      | \
// Eta, |  \
//      |   \
//      |    \
//     3|_____\1
//         Xi


// Boundary Condition Flag:
// 0 means There is no BC at this node
// 1 means Essential
// 2 means Natural

// Velocity BCs
// Essential:	Entrance, Top and Bottom
// Natural:	Exit

// Polymer
// Essential:	Entrance
// Natural:	Top, Bottom, Exit

#include "VelBioMesh2d.h"

void VelBioMesh2d(double dx, 
		  double dy, 
		  double XL,
		  double YB,
		  int NX, 
		  int NY,
		  int Nem, 
		  int & TotalNumNod, 
		  int & Npe, 
		  E_SDM & Glxy, 
		  E_ISDM & Nod,
		  E_ISDM & Nod_BC_Hor,
		  E_ISDM & Nod_BC_Ver,
		  E_ISDM & Bio_Nod_BC) {
    // Variable declarations
  int Eln1, Eln2;
  
  // The number of elements is double the number of intervals in the X direction because we are using trianglular elements
  NX = 2 * NX;
  NY = 2 * NY;
  
  // The step size between nodes is half because we are using quadratice interpolating polynomials
  dx = dx / 2;
  dy = dy / 2;
  
  // size of various structures
  // Number of nodes in the mesh
  TotalNumNod = (NY + 1) * (NX + 1);
  
  // Number of nodes in each element
  Npe = 6;
  
  // GLobal coordinates for each node point, row acccessed by its global node number and
  // column access is as follows, 0 is for x and 1 is for y
  Glxy.Shape(TotalNumNod,2);

  // Allocate the apprpriate amount of memory for each of these structure
  // There are 6 nodes per element and any one of them could be a BC of some kind
  Nod.Shape(Nem,Npe);
  Nod_BC_Hor.Shape(Nem,Npe);
  Nod_BC_Ver.Shape(Nem,Npe);
  Bio_Nod_BC.Shape(Nem,Npe);

  // global coordinate of nodes
  for (int i = 0; i <= TotalNumNod-1; i++)
  {
    // x-value, modded global node number by number of nodes per row, so the remainder gives the x position
    Glxy(i,0) = (i % (NX + 1)) * dx + XL;
    
    // y-value, the floor gives how many rows have been completed, which is the height of the node.
    Glxy(i,1) = floor(i / (NX + 1)) * dy + YB;
  }

  // connectivity matrix, global node numbers and nodes start at 1 not 0.
  // The number of elements per row is NX and the number of nodes per row is NX + 1
  // Each row of elements spans 3 rows of nodes but the top row of nodes is the bottom 
  // row of nodes for the next row of elements
  for (int j = 0; j <= NY/2-1; j++) 	// indexes through the rows
  {
    for (int i = 0; i <= NX/2-1; i++)	// indexes through the columns
    {
      Eln1 = j*NX + 2*(i + 1) - 1;	// lower element number,
      Eln2 = Eln1 + 1;			// upper element number

      // Refer to the diagram for an explanation of the physical structures that are being taken 
      // advantage of to make these assignments
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
      
      // If the upper element is at the end of a row
      // The velocity and polymer have natural BCs at the exit
      // We skip the last node on the first row of nodes.  This node gets taken care of by the
      // "bottom" conditional statement.  Any assigned NBCs on the corners get over written as we progress down
      // through the conditional statements where appropriate.
      if( Eln2 % NX == 0)
      {
// 	if(Eln1 == (NX - 1)) // If the lower element is at the end of the bottom row
// 	{
// 	  //Lower Element
// 	  // Only nodes 2 and 5 so that we dont assign natural BCs to the lower corner
// 	  // because it is an essential BC node.
// 	  Nod_BC_Hor(Eln1-1,2) = 2; // Flag = 2 means Natural BC
// 	  Nod_BC_Hor(Eln1-1,4) = 2;
// 	  
// 	  Bio_Nod_BC(Eln1-1,2) = 2; // Flag = 2 means Natural BC
// 	  Bio_Nod_BC(Eln1-1,4) = 2;
// 	  
// 	  //Upper Element
// 	  Nod_BC_Hor(Eln2-1,1) = 2;
// 	  
// 	  Bio_Nod_BC(Eln2-1,1) = 2;
// 	  
// 	  //Lower Element
// 	  Nod_BC_Ver(Eln1-1,2) = 1; // Flag = 2 means Natural BC
// 	  Nod_BC_Ver(Eln1-1,4) = 1;
// 	  
// 	  //Upper Element
// 	  Nod_BC_Ver(Eln2-1,1) = 1;
// 	}
// 	else if(Eln1 == (NX * (NY / 2) - 1)) // If the lower element is at the end of the top row
// 	{
// 	  // The upper and lower elements don't have a Natural BCs at the end of the top row since
// 	  // nodes 2 and 3 respectively are essential BCs.
// 	  //Lower Element
// 	  Nod_BC_Hor(Eln1-1,4) = 2;
// 	  
// 	  Bio_Nod_BC(Eln1-1,4) = 2;
// 	  
// 	  //Lower Element
// 	  Nod_BC_Ver(Eln1-1,4) = 1;
// 	}
// 	else
// 	{
	  //Lower Element
	  // Nod_BC_Hor(Eln1-1,1) = 2;
	  Nod_BC_Hor(Eln1-1,2) = 2; // Flag = 2 means Natural BC
	  Nod_BC_Hor(Eln1-1,4) = 2;
	  
	  // Bio_Nod_BC(Eln1-1,1) = 2;
	  Bio_Nod_BC(Eln1-1,2) = 2; // Flag = 2 means Natural BC
	  Bio_Nod_BC(Eln1-1,4) = 2;
	  
	  //Upper Element
	  Nod_BC_Hor(Eln2-1,1) = 2;
	  
	  Bio_Nod_BC(Eln2-1,1) = 2;
	  
	  //Lower Element
	  // Nod_BC_Ver(Eln1-1,1) = 1;
	  Nod_BC_Ver(Eln1-1,2) = 1; // Flag = 2 means Natural BC
	  Nod_BC_Ver(Eln1-1,4) = 1;
	  
	  //Upper Element
	  Nod_BC_Ver(Eln2-1,1) = 1;
// 	}
      }
      
      // Boundary Condition assignment, specific to the mesh I generated
      // If the element is on the bottom row
      if( Eln2 <= NX ) // ( Nod(Eln1-1,0) >= 1) & (Nod(Eln1-1,0) <= NX) )
      {
	//Lower Element
	// The no slip lower boundary condition mean that the velocity
	// is specified here.
	Nod_BC_Hor(Eln1-1,0) = 1; // Flag = 1 means Essential BC
	Nod_BC_Hor(Eln1-1,1) = 1;
	Nod_BC_Hor(Eln1-1,3) = 1;
	
	// The polymer cannot flow through the lower boundary.
	Bio_Nod_BC(Eln1-1,0) = 2; // Flag = 2 means Natural BC
	Bio_Nod_BC(Eln1-1,1) = 2;
	Bio_Nod_BC(Eln1-1,3) = 2;
	
	//Upper Element
	Nod_BC_Hor(Eln2-1,0) = 1;
	
	Bio_Nod_BC(Eln2-1,0) = 2;
	
	//Lower Element
	Nod_BC_Ver(Eln1-1,0) = 1;
	Nod_BC_Ver(Eln1-1,1) = 1;
	Nod_BC_Ver(Eln1-1,3) = 1;
	
	//Upper Element
	Nod_BC_Ver(Eln2-1,0) = 1;
      }
      
      // Nod 3 from the lower element On the Top row
      if(Eln1 > Nem - NX) // ( Nod(Eln1-1,2) >= (NX + 1) * NY ) & ( Nod(Eln1-1,2) <= TotalNumNod ) )
      {
	//Lower Element
	// The no slip upper boundary condition mean that the velocity
	// is specified here.
	Nod_BC_Hor(Eln1-1,2) = 1; // Flag = 1 means Essential BC
	
	// The polymer cannot flow through the upper boundary.
	Bio_Nod_BC(Eln1-1,2) = 2; // Flag = 2 means Natural BC
	
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
      // By doing the front after the top and bottom we ensure that the polymer BC on the bottom and top front corners
      // are labeled as essential BCs
      if( Eln1 % NX == 1 ) // (Nod(Eln1-1,0)) % (NX + 1) == 1 )
      {
	//Lower Element
	// The inflow for the velocity is specified.
	Nod_BC_Hor(Eln1-1,0) = 1; // Flag = 1 means Essential BC
	
	// The inflow for the polymer is specified.
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
    }
  }
}