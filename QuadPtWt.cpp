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
  else if (NTRIQUAD == 12)
  {
    Tri_Quad_Pt.Shape(12,2);
    Tri_Quad_Wt.Shape(12,1);
    
    Tri_Quad_Pt(0,0) = 0.873821971016996;
    Tri_Quad_Pt(0,1) = 0.063089014491502;
    
    Tri_Quad_Pt(1,0) = 0.063089014491502;
    Tri_Quad_Pt(1,1) = 0.873821971016996;
    
    Tri_Quad_Pt(2,0) = 0.063089014491502;
    Tri_Quad_Pt(2,1) = 0.063089014491502;
    
    Tri_Quad_Pt(3,0) = 0.501426509658179;
    Tri_Quad_Pt(3,1) = 0.249286745170910;
    
    Tri_Quad_Pt(4,0) = 0.501426509658179;
    Tri_Quad_Pt(4,1) = 0.249286745170910;
    
    Tri_Quad_Pt(5,0) = 0.249286745170910;
    Tri_Quad_Pt(5,1) = 0.249286745170910;
    
    Tri_Quad_Pt(6,0) = 0.636502499121399;
    Tri_Quad_Pt(6,1) = 0.310352451033785;
    
    Tri_Quad_Pt(7,0) = 0.310352451033785;
    Tri_Quad_Pt(7,1) = 0.636502499121399;
    
    Tri_Quad_Pt(8,0) = 0.636502499121399;
    Tri_Quad_Pt(8,1) = 0.053145049844816;
    
    Tri_Quad_Pt(9,0) = 0.053145049844816;
    Tri_Quad_Pt(9,1) = 0.636502499121399;
    
    Tri_Quad_Pt(10,0) = 0.310352451033785;
    Tri_Quad_Pt(10,1) = 0.053145049844816;
    
    Tri_Quad_Pt(11,0) = 0.053145049844816;
    Tri_Quad_Pt(11,1) = 0.310352451033785;
    
    Tri_Quad_Wt(0,0) = 0.050844906370207;
    Tri_Quad_Wt(1,0) = 0.050844906370207;
    Tri_Quad_Wt(2,0) = 0.050844906370207;
    Tri_Quad_Wt(3,0) = 0.116786275726379;
    Tri_Quad_Wt(4,0) = 0.116786275726379;
    Tri_Quad_Wt(5,0) = 0.116786275726379;
    Tri_Quad_Wt(6,0) = 0.082851075618374;
    Tri_Quad_Wt(7,0) = 0.082851075618374;
    Tri_Quad_Wt(8,0) = 0.082851075618374;
    Tri_Quad_Wt(9,0) = 0.082851075618374;
    Tri_Quad_Wt(10,0) = 0.082851075618374;
    Tri_Quad_Wt(11,0) = 0.082851075618374;
  }
  else if (NTRIQUAD == 13)
  {
    Tri_Quad_Pt.Shape(13,2);
    Tri_Quad_Wt.Shape(13,1);
    
    Tri_Quad_Pt(0,0) = 1.0/3.0;
    Tri_Quad_Pt(0,1) = 1.0/3.0;
    
    Tri_Quad_Pt(1,0) = 0.479308067841923;
    Tri_Quad_Pt(1,1) = 0.260345966079038;
    
    Tri_Quad_Pt(2,0) = 0.260345966079038;
    Tri_Quad_Pt(2,1) = 0.479308067841923;
    
    Tri_Quad_Pt(3,0) = 0.260345966079038;
    Tri_Quad_Pt(3,1) = 0.260345966079038;
    
    Tri_Quad_Pt(4,0) = 0.869739794195568;
    Tri_Quad_Pt(4,1) = 0.065130102902216;
    
    Tri_Quad_Pt(5,0) = 0.065130102902216;
    Tri_Quad_Pt(5,1) = 0.869739794195568;
    
    Tri_Quad_Pt(6,0) = 0.065130102902216;
    Tri_Quad_Pt(6,1) = 0.065130102902216;
    
    Tri_Quad_Pt(7,0) = 0.638444188569809;
    Tri_Quad_Pt(7,1) = 0.312865496004875;
    
    Tri_Quad_Pt(8,0) = 0.312865496004875;
    Tri_Quad_Pt(8,1) = 0.638444188569809;
    
    Tri_Quad_Pt(9,0) = 0.638444188569809;
    Tri_Quad_Pt(9,1) = 0.048690315425316;
    
    Tri_Quad_Pt(10,0) = 0.048690315425316;
    Tri_Quad_Pt(10,1) = 0.638444188569809;
    
    Tri_Quad_Pt(11,0) = 0.048690315425316;
    Tri_Quad_Pt(11,1) = 0.312865496004875;
    
    Tri_Quad_Pt(12,0) = 0.312865496004875;
    Tri_Quad_Pt(12,1) = 0.048690315425316;
    
    Tri_Quad_Wt(0,0) = -0.149570044467670;
    Tri_Quad_Wt(1,0) = 0.175615257433204;
    Tri_Quad_Wt(2,0) = 0.175615257433204;
    Tri_Quad_Wt(3,0) = 0.175615257433204;
    Tri_Quad_Wt(4,0) = 0.053347235608839;
    Tri_Quad_Wt(5,0) = 0.053347235608839;
    Tri_Quad_Wt(6,0) = 0.053347235608839;
    Tri_Quad_Wt(7,0) = 0.077113760890257;
    Tri_Quad_Wt(8,0) = 0.077113760890257;
    Tri_Quad_Wt(9,0) = 0.077113760890257;
    Tri_Quad_Wt(10,0) = 0.077113760890257;
    Tri_Quad_Wt(11,0) = 0.077113760890257;
    Tri_Quad_Wt(12,0) = 0.077113760890257;
  }
}