// Requires NX equals NY

#include "InitialBio.h"

void InitialBio(int NX,
		E_SDM Glxy,
		double *Bio) 
{
  double x, y;
  int count;
  
  count = 0;
  
  for (int i = 0; i <= 2 * NX; i++) // Number nodes across a physical row and column are 2*NX+1
  {
    for(int j = 0; j <= 2 * NX; j++)
    {
      x = Glxy(count,0);
      y = Glxy(count,1);
      
      // Square of side length 1/2 centered on the x-axis.
//       if((x >= 1.0 / 4.0) & (x <= 3.0 / 4.0) & (y <= 1.0 / 2.0)) Bio(count) = 1.0;
      
      // Circle of radius 1/3 centered at (1/2,1/2 + 1/10).
//       if((x - 1.0 / 2.0) * (x - 1.0 / 2.0) + (y - 1.0 / 2.0 - 1.0 / 10.0) * (y - 1.0 / 2.0 - 1.0 / 10.0) < 1.0 / 9.0) Bio[count] = 1.0;
      
      // Circle of radius 1/2 centered at (1/2,1/2).
//       if((x - 1.0 / 2.0) * (x - 1.0 / 2.0) + (y - 1.0 / 2.0) * (y - 1.0 / 2.0) <= 1.0 / 4.0) Bio[count] = 1.0;
      
      // Half circle of radius 1/2 centered at (1/2,0).
//       if((x - 1.0 / 2.0) * (x - 1.0 / 2.0) + y * y <= 1.0 / 16.0) Bio[count] = 0.1;
      
      // Fat Cylinder
      if((x >= 1.0 / 8.0) && (x <= 3.0 / 8.0)  && (y <= sqrt(1.0 / 64.0 - (x - 1.0 / 4.0) * (x - 1.0 / 4.0)) + 1.0 / 8.0)) Bio[count] = 0.1;
      
      // Narrow tall cylinder to capture movement sooner
//       if((x >= 3.0 / 16.0) && (x <= 5.0 / 16.0)  && (y <= sqrt(1.0 / 256.0 - (x - 1.0 / 4.0) * (x - 1.0 / 4.0)) + 1.0 / 2.0)) Bio[count] = 0.1;
      
      // lollipop
      // Thin and tall, y less than 0.5
//       if((x >= 3.0 / 16.0) && (x <= 5.0 / 16.0)  && (y <= 0.5)) Bio[count];
      // Circle centered at (1/4,1/2) with radius 1/8
//       if((y <= sqrt(1.0 / 64.0 - (x - 1.0 / 4.0) * (x - 1.0 / 4.0)) + 1.0 / 4.0)) Bio[count];
      
      count++;
    }
  }
}