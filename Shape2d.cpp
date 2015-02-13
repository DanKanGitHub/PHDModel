// Shape function for 2D finite element method
// Input: Xi, Eta: natural coordinates
// 	  Npe: # of nodes per element
//        Elxy: coordinates of local nodes, array of size Npex2
//        FLAG: 1 = linear triangle,  2 = quadratic triangle
// Output: Sf: value of shape functions, size Npex1
//         Gdsf: Derivate w.r.t. global coordinates, size 2xNpe
//         DetJ: determinant of Jacobian

#include "Shape2d.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_SerialDenseVector E_SDV;

void Shape2d(double Xi,
	     double Eta,
	     E_SDM Elxy,
	     int Npe,
	     int FLAG,
	     E_SDV & Sf,
	     E_SDM & Gdsf,
	     double & DetJ)
{
  
  // code common to all types of elements
  E_SDM Dsf(2,Npe);  	// derivatives w.r.t. natural coordinates
  E_SDM Gj(2,2);	// Integral change of variables jacobian
  E_SDM Inv_Gj(2,2);	// We must convert back to our (x,y) coordinate space from barycentric coordinates

  for(int i = 0; i <= 1; i++)
  {
    for(int j = 0; j <= Npe - 1; j++)
    {
	Gdsf(i,j) = 0.0;
    }
    
    for(int j = 0; j <= 1; j++)
    {
	Gj(i,j) = 0.0;
    }
  }
  
  if (FLAG == 1) // linear triangle element
  {
    Sf(0) = Xi;
    Sf(1) = Eta;
    Sf(2) = 1.0 - Xi - Eta;
    
    Dsf(0,0) = 1.0;
    Dsf(0,1) = 0.0;
    Dsf(0,2) = -1.0;
    Dsf(1,0) = 0.0;
    Dsf(1,1) = 1.0;
    Dsf(1,2) = -1.0;
  }
  else if (FLAG == 2) // Quadratic triangular element
  {
    // pg 531 L1 = Xi, L2 = Eta, L3 = 1-Xi-Eta.  Book has wrong order?
    Sf(0) = 2.0 * Xi * Xi - Xi;
    Sf(3) = 4.0 * Xi * Eta;
    Sf(1) = 2.0 * Eta * Eta - Eta;
    Sf(4) = 4.0 * (Eta - Xi * Eta - Eta * Eta);
    Sf(2) = 2.0 * (1.0 - Xi - Eta) * (1.0 - Xi - Eta) - (1.0 - Xi - Eta);
    Sf(5) = 4.0 * (Xi - Xi * Xi - Eta * Xi);

    Dsf(0,0) = 4.0 * Xi - 1.0;
    Dsf(0,3) = 4.0 * Eta;
    Dsf(0,1) = 0.0;
    Dsf(0,4) = -4.0 * Eta;
    Dsf(0,2) = -4.0 * (1.0 - Xi - Eta) + 1.0;
    Dsf(0,5) = 4.0 * (1.0 - 2.0 * Xi - Eta);

    Dsf(1,0) = 0.0;
    Dsf(1,3) = 4.0 * Xi;
    Dsf(1,1) = 4.0 * Eta - 1.0;
    Dsf(1,4) = 4.0 * (1.0 - Xi - 2.0 * Eta);
    Dsf(1,2) = -4.0 * (1.0 - Xi - Eta) + 1.0;
    Dsf(1,5) = -4.0 * Xi;
  }

  // Jacobian  
  for(int j = 0; j <= 1; j++)
  {
    for(int i = 0; i <= 1; i++)
    {
      for(int k = 0; k <= Npe - 1; k++)
      {
	Gj(i,j) += Dsf(i,k) * Elxy(k,j);
      }
    }
  }

  DetJ = Gj(0,0)*Gj(1,1) - Gj(0,1)*Gj(1,0);

  Inv_Gj(0,0) = 1.0 / DetJ * Gj(1,1);
  Inv_Gj(0,1) = -1.0 / DetJ * Gj(0,1);
  Inv_Gj(1,0) = -1.0 / DetJ * Gj(1,0);
  Inv_Gj(1,1) = 1.0 / DetJ * Gj(0,0);

  for(int j = 0; j <= Npe-1; j++)
  {
    for(int i = 0; i <= 1; i++)
    {
      for(int k = 0; k <= 1; k++)
      {
	Gdsf(i,j) += Inv_Gj(i,k) * Dsf(k,j);
      }
    }
  }
}
