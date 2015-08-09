// Departure Feet algorithm
// global coordinates in and out (that is, not barycentric coordinates)

#include "DepartureFoot.h"

void DepartureFoot(int Vel_Flag, 
		   int Nem, 
		   int Vel_Nnm,
		   int NX,
		   double TIMESTEP,
		   double dx,
		   double dy,
		   double *Vel, 
		   E_SDM Vel_Glxy, 
		   int Vel_Npe, 
		   E_ISDM Vel_Nod, 
		   E_ISDM Ele_Neigh,
		   E_SDV Ini_Foot,
		   E_SDV & Depart_Foot,
		   int & New_Ele) 
{
  
  // Velocity terms
  E_SDM El_Vel, Grad_Vel, Det_Inv_Grad_Vel, Vel_Elxy; // , Ini_Depart_Vel;
  E_SDV It_Vel; // iterative velocity
  double Det_Grad_Vel, Trace_Grad_Vel;

  // Other
  E_SDM Gdsf, TempMat;
  E_SDV Y_Norm, X, Y_New, Y_Old, Y_New_Xi_Eta, Sf, TempVec;
  int Inod, Cur_Ele, Old_Ele, Ele_Flag;
  double Xi, Eta, DetJ, Coeff, Feet_Tol, TOL = 0.00001; // Same tolerance as in main
  double alpha1, alpha2, alpha3, x, y, x1, x2, x3, y1, y2, y3, Two_Area;
  double xDbl, yDbl;
  
  // Initialize
  // Velocity terms
  El_Vel.Shape(2,Vel_Npe);
  It_Vel.Size(2);
  Det_Inv_Grad_Vel.Shape(2,2);
  Grad_Vel.Shape(2,2);
  Vel_Elxy.Shape(Vel_Npe,2);

  // Other terms
  X.Size(2);
  Y_New.Size(2);
  Y_Old.Size(2);
  Y_New_Xi_Eta.Size(2);
  Y_Norm.Size(2);
  Sf.Size(Vel_Npe);  		// value of shape functions at (xi,eta)
  Gdsf.Shape(2,Vel_Npe); 	// derivatives w.r.t. global cooridinates
  TempMat.Shape(2,2);
  TempVec.Size(2);

  // Feet Search Variables
  int xInt, yInt;
  double slope; //xRemainder, yRemainder, 
  
  // Initialize the while loop
  // Element velocity info
  for (int i = 0; i <= Vel_Npe-1; i++) 	// get global coordinates of local nodes of element NE
  {
    Inod = Vel_Nod(Nem, i)-1;		// Global node number (minus one for indexing) of local node.
    Vel_Elxy(i,0) = Vel_Glxy(Inod, 0);	// x-coordinate
    Vel_Elxy(i,1) = Vel_Glxy(Inod, 1);	// y-coordinate
  
    El_Vel(0,i) = Vel[Inod];
    El_Vel(1,i) = Vel[Inod + Vel_Nnm];
  }

  alpha1 = Vel_Elxy(1,0) * Vel_Elxy(2,1) - Vel_Elxy(2,0) * Vel_Elxy(1,1); // x2 * y3 - x3 * y2
  alpha2 = Vel_Elxy(2,0) * Vel_Elxy(0,1) - Vel_Elxy(0,0) * Vel_Elxy(2,1); // x3 * y1 - x1 * y3
  alpha3 = Vel_Elxy(0,0) * Vel_Elxy(1,1) - Vel_Elxy(1,0) * Vel_Elxy(0,1); // x1 * y2 - x2 * y1
  
  Two_Area = alpha1 + alpha2 + alpha3;
  
  //Initial position and constant throughout this subroutine
  x = Ini_Foot(0);
  y = Ini_Foot(1);
  
  x1 = Vel_Elxy(0,0);
  x2 = Vel_Elxy(1,0);
  x3 = Vel_Elxy(2,0);
  y1 = Vel_Elxy(0,1);
  y2 = Vel_Elxy(1,1);
  y3 = Vel_Elxy(2,1);
  
  //Convert the initial position into barycentric coordinates
  Xi  = 1.0 / Two_Area * ((x - x3) * (y2 - y3) - (y - y3) * (x2 - x3));
  Eta = 1.0 / Two_Area * ((x - x1) * (y3 - y1) + (y - y1) * (x1 - x3));

  // X is the starting position vector, in barycentric coordinates, and is constant throughout the loop
  X(0) = Xi;
  X(1) = Eta;

  // Y_Old is the starting position vector, in cartesian coordinates.
  Y_Old(0) = x; // Initial Guess is at X
  Y_Old(1) = y;

  Shape2d(Xi,
	  Eta,
	  Vel_Elxy,
	  Vel_Npe,
	  Vel_Flag,
	  Sf,	// output
	  Gdsf,	// output
	  DetJ);	// output

  // Zero out data structures for use of += later.
  for(int i = 0; i <= 1; i++)
  {
    for(int j = 0; j <= 1; j++)
    {
      Grad_Vel(i,j) = 0.0;
    }
    It_Vel(i) = 0.0;
  }
  
  // Grad_Vel = El_Vel * Gdsf' (transpose)
  for(int i = 0; i <= 1; i++)
  {
    for(int j = 0; j <= 1; j++)
    {
      for(int k = 0; k <= Vel_Npe-1; k++)
      {
	Grad_Vel(i,j) += El_Vel(i,k) * Gdsf(j,k); // Gdsf transpose
      }
    }
  }
  
  // Iteratived velocity It_Vel = El_Vel * Sf
  for(int j = 0; j <= Vel_Npe - 1; j++)
  {
    It_Vel(0) += El_Vel(0,j) * Sf(j);
    It_Vel(1) += El_Vel(1,j) * Sf(j);
  }

//   Y_Old(0) = x - TIMESTEP * It_Vel(0); // Initial Guess is at X
//   Y_Old(1) = y - TIMESTEP * It_Vel(1);
  
  Det_Grad_Vel = Grad_Vel(0,0)*Grad_Vel(1,1) - Grad_Vel(0,1)*Grad_Vel(1,0);

  Trace_Grad_Vel = Grad_Vel(0,0) + Grad_Vel(1,1);

  // det(Inv_Grad_Vel)
  Coeff = 1.0 + TIMESTEP / 2.0 * Trace_Grad_Vel + TIMESTEP * TIMESTEP / 4.0 * Det_Grad_Vel;

  // The invere of gradu times it's determinant
  Det_Inv_Grad_Vel(0,0) = Grad_Vel(1,1);
  Det_Inv_Grad_Vel(0,1) = -Grad_Vel(0,1);
  Det_Inv_Grad_Vel(1,0) = -Grad_Vel(1,0);
  Det_Inv_Grad_Vel(1,1) = Grad_Vel(0,0);

  // Computing I - k/2 * det(gradu)*(gradu Inv)
  TempMat(0,0) = 1.0 + TIMESTEP / 2.0 * Det_Inv_Grad_Vel(0,0);
  TempMat(0,1) = -TIMESTEP / 2.0 * Det_Inv_Grad_Vel(0,1);
  TempMat(1,0) = -TIMESTEP / 2.0 * Det_Inv_Grad_Vel(1,0);
  TempMat(1,1) = 1.0 + TIMESTEP / 2.0 * Det_Inv_Grad_Vel(1,1);

  // computing yn + k * u((yn + x)/2) - x
  TempVec(0) = Y_Old(0) + TIMESTEP * It_Vel(0) - x;
  TempVec(1) = Y_Old(1) + TIMESTEP * It_Vel(1) - y;

  Y_New(0) = Y_Old(0) - 1.0 / Coeff * (TempMat(0,0) * TempVec(0) + TempMat(0,1) * TempVec(1));
  Y_New(1) = Y_Old(1) - 1.0 / Coeff * (TempMat(1,0) * TempVec(0) + TempMat(1,1) * TempVec(1));
  
  // Returns the departure foot to the domain
  if(Y_New(0) < 0)
  {
    
//     std::cout << "Y_New(0) = " << Y_New(0) << endl;
    
    Y_New(0) = 0;
  }
  
  if(Y_New(1) < 0)
  {
    Y_New(1) = 0;
  }
  
  if(Y_New(0) > 1)
  {
    Y_New(0) = 1;
  }
  
  if(Y_New(1) > 1)
  {
    Y_New(1) = 1;
  }
  
  // Compute the vector for determining the two norm of the differenece in iteraions
  Y_Norm(0) = Y_New(0) - Y_Old(0);
  Y_Norm(1) = Y_New(1) - Y_Old(1);

  // Compute the two norm
  Feet_Tol = Y_Norm.Norm2();// sqrt(Y_Norm(0,0)*Y_Norm(0,0) + Y_Norm(0,1)*Y_Norm(0,1))

  Cur_Ele = Nem + 1; // Element numbering starts at 1 but Ne starts at 0
  
  Old_Ele = Cur_Ele;
  New_Ele = -1;
  Ele_Flag = 1;

//       int Count = 1;
  
//   std::cout << "Before Depart_Foot while loop" << endl;
  
  int count = 0;
  int QWERTY;
  
  // Do not enter the while loop if already converged
  while (Feet_Tol > TOL)// && counter <= 10) // converges in at most 5 steps
  {
    
    count++;
    
    if(count >= 10) {
      
//       std::cout << "Cur_Ele = " << Cur_Ele << endl;
//       std::cout << "New_Ele = " << New_Ele << endl;
//       std::cout << "Old_Ele = " << Old_Ele << endl;
//       std::cout << "Count = " << count << endl;
//       
//       std::cout << "Y_Old(0) = " << Y_Old(0) << endl;
//       std::cout << "Y_Old(1) = " << Y_Old(1) << endl;
//       std::cout << "Y_New(0) = " << Y_New(0) << endl;
//       std::cout << "Y_New(1) = " << Y_New(1) << endl;
//       std::cout << "xInt = " << xInt << endl;
//       std::cout << "yInt = " << yInt << endl;
      
//       std::cin >> QWERTY;
      
      // Computer Stress at original point
      
      // Compute Stress at Old point
      
      // Compute Stress at new point
      
      // Determine the point with the Stress closest to the original Stress
      
      // Assign the associated element and point
      
      // break out of the while loop
      
      // Average Y_Old and Y_New
//       Y_New(0) = (Y_New(0) + Y_Old(0)) / 2;
//       Y_New(1) = (Y_New(1) + Y_Old(1)) / 2;
      
      count = 0;
      
      TIMESTEP = TIMESTEP / 10;
      
    }

    //Convert the new position to barycentric coordinates.
//     Y_New_Xi_Eta(0) = 1.0 / Two_Area * ((Y_New(0) - Vel_Elxy(2,0)) * (Vel_Elxy(1,1) - Vel_Elxy(2,1)) - 
// 		    (Y_New(1) - Vel_Elxy(2,1)) * (Vel_Elxy(1,0) - Vel_Elxy(2,0)));
//   
//     Y_New_Xi_Eta(1) = 1.0 / Two_Area * ((Y_New(0) - Vel_Elxy(0,0)) * (Vel_Elxy(2,1) - Vel_Elxy(0,1)) + 
// 		    (Y_New(1) - Vel_Elxy(0,1)) * (Vel_Elxy(0,0) - Vel_Elxy(2,0)));

//	// Integer division to get the row and column
    xInt = Y_New(0) / dx;
    yInt = Y_New(1) / dy;
    
    xDbl = floor(Y_New(0) / dx);
    yDbl = floor(Y_New(1) / dy);
    
//     std::cout << "Y_Old(0) = " << Y_Old(0) << endl;
//     std::cout << "Y_Old(1) = " << Y_Old(1) << endl;
//     std::cout << "Y_New(0) = " << Y_New(0) << endl;
//     std::cout << "Y_New(1) = " << Y_New(1) << endl;
//     std::cout << "xInt = " << xInt << endl;
//     std::cout << "yInt = " << yInt << endl;
//     std::cout << "Cur_Ele = " << Cur_Ele << endl;
    

//     xRemainder = Y_New(0) / dx;
//     yRemainder = Y_New(1) / dy;
    
    // I don't need this
//     slope = Y_New(1) / Y_New(0);
    
//     if ((xInt != NX) & (yInt != NX)) {
//       
//       // If xInt = yInt then compare to the line y = x.  If xInt + 1 = yInt then compare to the line y = x + deltax, etc.
//       // In general, let n = yInt - xInt and compare to the line y = x + n * deltax, where deltax = 1 / NX.
//       // Assigns the point to the lower element if it is on the diagonal line.
//       // If y is above the line
//       if (Y_New(1) > Y_New(0) + (yInt - xInt) / NX) {
// 	
// 	New_Ele = yInt * 2 * NX + 2 * xInt + 2; // Top element
// 	
//       } else {
// 	
// 	New_Ele = yInt * 2 * NX + 2 * xInt + 1; // bottom element
// 	
//       };
//     } else if ((xInt == NX) & (yInt != NX)) {
// 	
//       New_Ele = (yInt + 1) * 2 * NX - 1; // Always the lower element when the point is on the right boundary
// 	
//     } else if ((xInt != NX) & (yInt == NX)) {
// 
//       New_Ele = 2 * NX * (NX - 1) + 2 * (xInt + 1); // always the upper element when the point is on the upper boundary
// 
//     } else if ((xInt == NX) & (yInt == NX)) {
// 
//       New_Ele = 2 * NX * NX - 1; // lower element when the point is the upper right corner.
// 
//     };

    // if( abs(  ) < 0.000001)
    // if( Y_New(0) / dx != floor(Y_New(0) / dx) )
    
//     if ((fmod(Y_New(0), dx) != 0.0) & (fmod(Y_New(1), dy) != 0.0)) {
    if( (Y_New(0) / dx != floor(Y_New(0) / dx)) & (Y_New(1) / dx != floor(Y_New(1) / dy)) ) {
      
      // If xInt = yInt then compare to the line y = x.  If xInt + 1 = yInt then compare to the line y = x + deltax, etc.
      // In general, let n = yInt - xInt and compare to the line y = x + n * deltax, where deltax = 1 / NX.
      // Assigns the point to the lower element if it is on the diagonal line.
      // If y is above the line
      if (Y_New(1) > Y_New(0) + (yDbl - xDbl) * dx) {
	
	New_Ele = yInt * 2 * NX + 2 * xInt + 2; // Top element
	
      } else {
	
	New_Ele = yInt * 2 * NX + 2 * xInt + 1; // bottom element
	
      };
//     } else if ((fmod(Y_New(0), dx) == 0.0) & (fmod(Y_New(1), dy) != 0.0)) {
    } else if ( (Y_New(0) / dx == floor(Y_New(0) / dx)) & (Y_New(1) / dx != floor(Y_New(1) / dy)) ) {
      
      if(Y_New(0) != 0.0) {
	
	New_Ele = yInt * 2 * NX + 2 * xInt - 1; // Always the lower element when the point is on the right boundary
	
      } else {
	
	New_Ele = yInt * 2 * NX + 2; // On the front of the domain
	
      };
	
//     } else if ((fmod(Y_New(0), dx) != 0.0) & (fmod(Y_New(1), dy) == 0.0)) {
    } else if ((Y_New(0) / dx != floor(Y_New(0) / dx)) & (Y_New(1) / dx == floor(Y_New(1) / dy))) {
      
      if(Y_New(1) != 0.0) {

	New_Ele = 2 * NX * (yInt - 1) + 2 * (xInt + 1); // always the upper element when the point is on the upper boundary
	
      } else {
	
	New_Ele = 2 * xInt + 1; // On the bottom of the domain.
	
      };

//     } else if ((fmod(Y_New(0), dx) == 0.0) & (fmod(Y_New(1), dy) == 0.0)) {
    } else if ((Y_New(0) / dx == floor(Y_New(0) / dx)) & (Y_New(1) / dx == floor(Y_New(1) / dy))) {
      
      if((Y_New(0) != 0.0) & (Y_New(1) != 0.0)) {

	New_Ele = 2 * NX * (yInt - 1) + 2 * xInt - 1; // lower element when the point is the upper right corner.
	
      } else if((Y_New(0) == 0.0) & (Y_New(1) != 0.0)){
	
	New_Ele = 2 * NX * (yInt - 1) + 2;
      
      } else if((Y_New(0) != 0.0) & (Y_New(1) == 0.0)){
	
	New_Ele = 2 * xInt - 1;
      
      } else {

	New_Ele = 1; // The lower left corner.
	
      };

    };
    
//     std::cout << "New_Ele = " << New_Ele << endl;
    
//     int QWERTY;
//     std::cout << "Cur_Ele = " << Cur_Ele << endl;
    
//     while(Ele_Flag == 1) {
      // Checks which element the point is in using barycentric coordinates.
//       New_Ele = FeetSearch(Y_New_Xi_Eta,
// 			    Cur_Ele, // Old_Ele
// 			    Ele_Neigh);
//       
//       std::cout << "Before Departure Element" << endl;
//       std::cout << "New_Ele = " << New_Ele << endl;
//       std::cout << "Ele_Neigh = " << Ele_Neigh << endl;
//       std::cin >> QWERTY;
// 
//       if(Old_Ele == New_Ele) {
// 	Ele_Flag = 0;
//       } else {
// 	for (int i = 0; i <= Vel_Npe-1; i++) 	// get global coordinates of local nodes of element NE
// 	{
// 	  Inod = Vel_Nod(New_Ele - 1, i)-1;		// Global node number (minus one for indexing) of local node.
// 	  Vel_Elxy(i,0) = Vel_Glxy(Inod, 0);	// x-coordinate
// 	  Vel_Elxy(i,1) = Vel_Glxy(Inod, 1);	// y-coordinate
// 	}
// 	
// 	Y_New_Xi_Eta(0) = 1.0 / Two_Area * ((Y_New(0) - Vel_Elxy(2,0)) * (Vel_Elxy(1,1) - Vel_Elxy(2,1)) - 
// 		      (Y_New(1) - Vel_Elxy(2,1)) * (Vel_Elxy(1,0) - Vel_Elxy(2,0)));
//     
// 	Y_New_Xi_Eta(1) = 1.0 / Two_Area * ((Y_New(0) - Vel_Elxy(0,0)) * (Vel_Elxy(2,1) - Vel_Elxy(0,1)) + 
// 			(Y_New(1) - Vel_Elxy(0,1)) * (Vel_Elxy(0,0) - Vel_Elxy(2,0)));
// 	
// 	if(Y_New_Xi_Eta(0) < 0) {
// 	  
// 	  std::cout << "Y_New_Xi_Eta(0) = " << Y_New_Xi_Eta(0) << endl;
// 	  std::cout << "Y_New(0) = " << Y_New(0) << endl;
// 	  
// 	}
// 	
// 	Old_Ele =  New_Ele;
//       }; // end if
//     }; // end ehile

//     Ele_Flag = 1;
    
    if(New_Ele != Cur_Ele)
    {
      // If the foot changed elements then update velocity info.
      for (int i = 0; i <= Vel_Npe-1; i++) // get global coordinates of local nodes of element NE
      {
	Inod = Vel_Nod(New_Ele-1,i)-1;		// Global node number of local node.
	Vel_Elxy(i,0) = Vel_Glxy(Inod,0);	// x-coordinate of te new element
	Vel_Elxy(i,1) = Vel_Glxy(Inod,1);	// y-coordinate of the new element
  
	El_Vel(0,i) = Vel[Inod];
	El_Vel(1,i) = Vel[Inod + Vel_Nnm];
      }
      
      // Update element number
      Cur_Ele = New_Ele;
      
    }

    // Update position vector
    Y_Old = Y_New;

    // averages the initial position with the new position and converts to barycentric coordinates
    Xi = 1.0 / Two_Area * (((Y_New(0) + x) / 2.0 - Vel_Elxy(2,0)) * (Vel_Elxy(1,1) - Vel_Elxy(2,1)) - 
	    ((Y_New(1) + y) / 2.0 - Vel_Elxy(2,1)) * (Vel_Elxy(1,0) - Vel_Elxy(2,0)));
  
    Eta = 1.0 / Two_Area * (((Y_New(0) + x) / 2.0 - Vel_Elxy(0,0)) * (Vel_Elxy(2,1) - Vel_Elxy(0,1)) + 
	    ((Y_New(1) + y) / 2.0 - Vel_Elxy(0,1)) * (Vel_Elxy(0,0) - Vel_Elxy(2,0)));

    // Recompute the shape function info using the new element info
    Shape2d(Xi, 
	    Eta, 
	    Vel_Elxy, 
	    Vel_Npe, 
	    Vel_Flag, 
	    Sf, 	// output
	    Gdsf,	// output
	    DetJ);	// output

    for(int i = 0; i <= 1; i++)
    {
      for(int j = 0; j <= 1; j++)
      {
	Grad_Vel(i,j) = 0.0;
      }
      It_Vel(i) = 0.0;
    }
  
    for(int i = 0; i <= 1; i++)
    {
      for(int j = 0; j <= 1; j++)
      {
	for(int k = 0; k <= Vel_Npe-1; k++)
	{
	  Grad_Vel(i,j) += El_Vel(i,k) * Gdsf(j,k);
	}
      }
    }

    for(int j = 0; j <= Vel_Npe - 1; j++)
    {
      It_Vel(0) += El_Vel(0,j) * Sf(j);
      It_Vel(1) += El_Vel(1,j) * Sf(j);
    }
    
    Det_Grad_Vel = Grad_Vel(0,0)*Grad_Vel(1,1) - Grad_Vel(0,1)*Grad_Vel(1,0);
    
    Trace_Grad_Vel = Grad_Vel(0,0) + Grad_Vel(1,1);
  
    Coeff = 1.0 + TIMESTEP / 2.0 * Trace_Grad_Vel + TIMESTEP * TIMESTEP / 4.0 * Det_Grad_Vel;

    // The invere of gradu times it's determinant
    Det_Inv_Grad_Vel(0,0) = Grad_Vel(1,1);
    Det_Inv_Grad_Vel(0,1) = -Grad_Vel(0,1);
    Det_Inv_Grad_Vel(1,0) = -Grad_Vel(1,0);
    Det_Inv_Grad_Vel(1,1) = Grad_Vel(0,0);

    // yold = x for the initializing step I - k/2 * det(gradu)*(gradu Inv)
    TempMat(0,0) = 1.0 + TIMESTEP / 2.0 * Det_Inv_Grad_Vel(0,0);
    TempMat(0,1) = -TIMESTEP / 2.0 * Det_Inv_Grad_Vel(0,1);
    TempMat(1,0) = -TIMESTEP / 2.0 * Det_Inv_Grad_Vel(1,0);
    TempMat(1,1) = 1.0 + TIMESTEP / 2.0 * Det_Inv_Grad_Vel(1,1);

    // yn + k * u((yn + x)/2) - x
    TempVec(0) = Y_Old(0) + TIMESTEP * It_Vel(0) - x;
    TempVec(1) = Y_Old(1) + TIMESTEP * It_Vel(1) - y;

    Y_New(0) = Y_Old(0) - 1.0 / Coeff * (TempMat(0,0) * TempVec(0) + TempMat(0,1) * TempVec(1));
    Y_New(1) = Y_Old(1) - 1.0 / Coeff * (TempMat(1,0) * TempVec(0) + TempMat(1,1) * TempVec(1));

    // Keep the departure foot in the domain
    if(Y_New(0) < 0)
    {
      Y_New(0) = 0;
    }
  
    if(Y_New(1) < 0)
    {
      Y_New(1) = 0;
    }
  
    if(Y_New(0) > 1)
    {
      Y_New(0) = 1;
    }
  
    if(Y_New(1) > 1)
    {
      Y_New(1) = 1;
    }
    
    // Compute the vector for determining the two norm of the differenece in iteraions
    Y_Norm(0) = Y_New(0) - Y_Old(0);
    Y_Norm(1) = Y_New(1) - Y_Old(1);

    // Compute the two norm
    Feet_Tol = Y_Norm.Norm2();

  } // while
  
//   std::cout << "After Depart_Foot while loop" << endl;

  //Convert the converged departure foot to barycentric coordinates
//   Y_New_Xi_Eta(0) = 1.0 / Two_Area * ((Y_New(0) - Vel_Elxy(2,0)) * (Vel_Elxy(1,1) - Vel_Elxy(2,1)) - 
// 		    (Y_New(1) - Vel_Elxy(2,1)) * (Vel_Elxy(1,0) - Vel_Elxy(2,0)));
//   
//   Y_New_Xi_Eta(1) = 1.0 / Two_Area * ((Y_New(0) - Vel_Elxy(0,0)) * (Vel_Elxy(2,1) - Vel_Elxy(0,1)) + 
// 		    (Y_New(1) - Vel_Elxy(0,1)) * (Vel_Elxy(0,0) - Vel_Elxy(2,0)));
  
  xInt = Y_New(0) / dx;
  yInt = Y_New(1) / dy;
  
  xDbl = floor(Y_New(0) / dx);
  yDbl = floor(Y_New(1) / dy);
  
//   if ((xInt != NX) & (yInt != NX)) {
//     
    // If xInt = yInt then compare to the line y = x.  If xInt + 1 = yInt then compare to the line y = x + deltax, etc.
    // In general, let n = yInt - xInt and compare to the line y = x + n * deltax, where deltax = 1 / NX.
    // Assigns the point to the lower element if it is on the diagonal line.
    // If y is above the line
//     if (Y_New(1) > Y_New(0) + (yInt - xInt) / NX) {
//       
//       New_Ele = yInt * 2 * NX + 2 * xInt + 2; // Top element
//       
//     } else {
//       
//       New_Ele = yInt * 2 * NX + 2 * xInt + 1; // bottom element
//       
//     };
//   } else if ((xInt == NX) & (yInt != NX)) {
//       
//     New_Ele = (yInt + 1) * 2 * NX - 1; // Always the lower element when the point is on the right boundary
//       
//   } else if ((xInt != NX) & (yInt == NX)) {
// 
//     New_Ele = 2 * NX * (NX - 1) + 2 * (xInt + 1); // always the upper element when the point is on the upper boundary
// 
//   } else if ((xInt == NX) & (yInt == NX)) {
// 
//     New_Ele = 2 * NX * NX - 1; // lower element when the point is the upper right corner.
// 
//   };

  //     if ((fmod(Y_New(0), dx) != 0.0) & (fmod(Y_New(1), dy) != 0.0)) {
    if( (Y_New(0) / dx != floor(Y_New(0) / dx)) & (Y_New(1) / dx != floor(Y_New(1) / dy)) ) {
      
      // If xInt = yInt then compare to the line y = x.  If xInt + 1 = yInt then compare to the line y = x + deltax, etc.
      // In general, let n = yInt - xInt and compare to the line y = x + n * deltax, where deltax = 1 / NX.
      // Assigns the point to the lower element if it is on the diagonal line.
      // If y is above the line
      if (Y_New(1) > Y_New(0) + (yDbl - xDbl) * dx) {
	
	New_Ele = yInt * 2 * NX + 2 * xInt + 2; // Top element
	
      } else {
	
	New_Ele = yInt * 2 * NX + 2 * xInt + 1; // bottom element
	
      };
//     } else if ((fmod(Y_New(0), dx) == 0.0) & (fmod(Y_New(1), dy) != 0.0)) {
    } else if ( (Y_New(0) / dx == floor(Y_New(0) / dx)) & (Y_New(1) / dx != floor(Y_New(1) / dy)) ) {
      
      if(Y_New(0) != 0.0) {
	
	New_Ele = yInt * 2 * NX + 2 * xInt - 1; // Always the lower element when the point is on the right boundary
	
      } else {
	
	New_Ele = yInt * 2 * NX + 2; // On the front of the domain
	
      };
	
//     } else if ((fmod(Y_New(0), dx) != 0.0) & (fmod(Y_New(1), dy) == 0.0)) {
    } else if ((Y_New(0) / dx != floor(Y_New(0) / dx)) & (Y_New(1) / dx == floor(Y_New(1) / dy))) {
      
      if(Y_New(1) != 0.0) {

	New_Ele = 2 * NX * (yInt - 1) + 2 * (xInt + 1); // always the upper element when the point is on the upper boundary
	
      } else {
	
	New_Ele = 2 * xInt + 1; // On the bottom of the domain.
	
      };

//     } else if ((fmod(Y_New(0), dx) == 0.0) & (fmod(Y_New(1), dy) == 0.0)) {
    } else if ((Y_New(0) / dx == floor(Y_New(0) / dx)) & (Y_New(1) / dx == floor(Y_New(1) / dy))) {
      
      if((Y_New(0) != 0.0) & (Y_New(1) != 0.0)) {

	New_Ele = 2 * NX * (yInt - 1) + 2 * xInt - 1; // lower element when the point is the upper right corner.
	
      } else if((Y_New(0) == 0.0) & (Y_New(1) != 0.0)){
	
	New_Ele = 2 * NX * (yInt - 1) + 2;
      
      } else if((Y_New(0) != 0.0) & (Y_New(1) == 0.0)){
	
	New_Ele = 2 * xInt - 1;
      
      } else {

	New_Ele = 1; // The lower left corner.
	
      };

    };

//   std::cout << "Before Departure Element" << endl;
//   std::cout << "Old_Ele = " << Old_Ele << endl;

//   while(Ele_Flag == 1) {
// 
    // Checks which element the point is in using barycentric coordinates.
//     New_Ele = FeetSearch(Y_New_Xi_Eta,
// 			  Cur_Ele, // Old_Ele
// 			  Ele_Neigh);
// 
//     if(Old_Ele == New_Ele) {
//       Ele_Flag = 0;
//     } else {
//       for (int i = 0; i <= Vel_Npe-1; i++) 	// get global coordinates of local nodes of element NE
//       {
// 	Inod = Vel_Nod(New_Ele - 1, i)-1;		// Global node number (minus one for indexing) of local node.
// 	Vel_Elxy(i,0) = Vel_Glxy(Inod, 0);	// x-coordinate
// 	Vel_Elxy(i,1) = Vel_Glxy(Inod, 1);	// y-coordinate
//       }
//       
//       Y_New_Xi_Eta(0) = 1.0 / Two_Area * ((Y_New(0) - Vel_Elxy(2,0)) * (Vel_Elxy(1,1) - Vel_Elxy(2,1)) - 
// 		    (Y_New(1) - Vel_Elxy(2,1)) * (Vel_Elxy(1,0) - Vel_Elxy(2,0)));
//   
//       Y_New_Xi_Eta(1) = 1.0 / Two_Area * ((Y_New(0) - Vel_Elxy(0,0)) * (Vel_Elxy(2,1) - Vel_Elxy(0,1)) + 
// 		      (Y_New(1) - Vel_Elxy(0,1)) * (Vel_Elxy(0,0) - Vel_Elxy(2,0)));
//       
//       Old_Ele =  New_Ele;
//     };
//     
//   };

  Depart_Foot(0) = Y_New(0);
  Depart_Foot(1) = Y_New(1);
  
//   std::cout << "Depart_Foot(0) = " << Depart_Foot(0) << endl;
//   std::cout << "Depart_Foot(1) = " << Depart_Foot(1) << endl;
//   std::cout << "New_Ele = " << New_Ele << endl;
//   
//   int QWERTY;
//   std::cin >> QWERTY;
  
}
