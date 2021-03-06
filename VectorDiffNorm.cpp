// Computes the 2-norm of the difference between Vector1 and Vector2

#include "VectorDiffNorm.h"

double VectorDiffNorm(double *Vector1, 
		      double *Vector2, 
		      int VectorLength)
{
  
  double y = 0.0;
  
//   cout << "VectorLength = " << VectorLength << endl;
  
  for(int i = 0; i < VectorLength; i++)
  {
    
    y += (Vector1[i] - Vector2[i]) * (Vector1[i] - Vector2[i]);
    
  }
  
//   cout << "y = " << y << endl;
  
//   y = sqrt(y);

  return sqrt(y); //y;
}