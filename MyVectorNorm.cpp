// Computes the 2-norm of a vector.

#include "MyVectorNorm.h"

double MyVectorNorm(double* Vector1, 
		  double* Vector2, 
		  int VectorLength)
{
  
  double y = 0.0;
  
  std::cout << "VectorLength = " << VectorLength << std::endl;
  
  for(int i = 0; i < VectorLength; i++)
  {
    
    y = y + (Vector1[i] - Vector2[i]) * (Vector1[i] - Vector2[i]);
    
  }
  
  std::cout << "y = " << y << endl;
  
  y = sqrt(y);

  return y;
}