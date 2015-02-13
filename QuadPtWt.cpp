// pg 559 FEM text
#include "QuadPtWt.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_SerialDenseVector E_SDV;

void QuadPtWt(int NTRIQUAD, 
	      E_SDM & Tri_Quad_Pt, 
	      E_SDV & Tri_Quad_Wt)
{
  
  if (NTRIQUAD == 1)
  {
    Tri_Quad_Pt.Shape(2,1);
    Tri_Quad_Wt.Shape(1,1);
   
    Tri_Quad_Pt(0,0) = 1.0/3.0;
    Tri_Quad_Pt(1,0) = 1.0/3.0;
    
    Tri_Quad_Wt(0,0) = 1.0;
  }
  else if (NTRIQUAD == 3)
  {
    Tri_Quad_Pt.Shape(3,2);
    Tri_Quad_Wt.Shape(3,1);
    
    Tri_Quad_Pt(0,0) = 1.0/2.0;
    Tri_Quad_Pt(0,1) = 0.0;
    Tri_Quad_Pt(1,0) = 1.0/2.0;
    Tri_Quad_Pt(1,1) = 1.0/2.0;
    Tri_Quad_Pt(2,0) = 0.0;
    Tri_Quad_Pt(2,1) = 1.0/2.0;
    
    Tri_Quad_Wt(0,0) = 1.0/3.0;
    Tri_Quad_Wt(1,0) = 1.0/3.0;
    Tri_Quad_Wt(2,0) = 1.0/3.0;
  }
  else if (NTRIQUAD == 4)
  {
    Tri_Quad_Pt.Shape(4,2);
    Tri_Quad_Wt.Shape(4,1);
    
    Tri_Quad_Pt(0,0) = 1.0/3.0;
    Tri_Quad_Pt(0,1) = 1.0/3.0;
    Tri_Quad_Pt(1,0) = 0.6;
    Tri_Quad_Pt(1,1) = 0.2;
    Tri_Quad_Pt(2,0) = 0.2;
    Tri_Quad_Pt(2,1) = 0.6;
    Tri_Quad_Pt(3,0) = 0.2;
    Tri_Quad_Pt(3,1) = 0.2;
    
    Tri_Quad_Wt(0,0) = -27.0/48.0;
    Tri_Quad_Wt(1,0) = 25.0/48.0;
    Tri_Quad_Wt(2,0) = 25.0/48.0;
    Tri_Quad_Wt(3,0) = 25.0/48.0;
  }
  else if (NTRIQUAD == 7)
  {
    Tri_Quad_Pt.Shape(7,2);
    Tri_Quad_Wt.Shape(7,1);
    
    Tri_Quad_Pt(0,0) = 1.0/3.0;
    Tri_Quad_Pt(0,1) = 1.0/3.0;
    Tri_Quad_Pt(1,0) = 0.797426985353;
    Tri_Quad_Pt(1,1) = 0.101286507323;
    Tri_Quad_Pt(2,0) = 0.101286507323;
    Tri_Quad_Pt(2,1) = 0.797426985353;
    Tri_Quad_Pt(3,0) = 0.101286507323;
    Tri_Quad_Pt(3,1) = 0.101286507323;
    Tri_Quad_Pt(4,0) = 0.059715871789;
    Tri_Quad_Pt(4,1) = 0.470142064105;
    Tri_Quad_Pt(5,0) = 0.470142064105;
    Tri_Quad_Pt(5,1) = 0.059715871789;
    Tri_Quad_Pt(6,0) = 0.470142064105;
    Tri_Quad_Pt(6,1) = 0.470142064105;
    
    
    Tri_Quad_Wt(0,0) = 0.225;
    Tri_Quad_Wt(1,0) = 0.125939180544;
    Tri_Quad_Wt(2,0) = 0.125939180544;
    Tri_Quad_Wt(3,0) = 0.125939180544;
    Tri_Quad_Wt(4,0) = 0.132394152788;
    Tri_Quad_Wt(5,0) = 0.132394152788;
    Tri_Quad_Wt(6,0) = 0.132394152788;
  }
}