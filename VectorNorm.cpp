// Computes the 2-norm of Vector1

#include "VectorNorm.h"

double VectorNorm(double *Vector1, 
		  int VectorLength)
{
  
  double y = 0.0;
  
//   cout << "VectorLength = " << VectorLength << endl;
  
  for(int i = 0; i < VectorLength; i++)
  {
    
    y += Vector1[i] * Vector1[i];
    
  }

  return sqrt(y); //y;
}