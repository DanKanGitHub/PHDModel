// Gaussian Filter Subroutine
// Takes in a matrix of data, pads the sides and runs the filter

// Matlab's Gaussian Filter sigma = 0.5, h = fspecial('gaussian', hsize, sigma)
//    0.011343736558495   0.083819505802211   0.011343736558495
//    0.083819505802211   0.619347030557177   0.083819505802211
//    0.011343736558495   0.083819505802211   0.011343736558495

void GaussianFilter(int NX, int NY, double & Bio_Data)
{
  E_SDM BioOutPut;
  
  for(int i = 0; i <= 2 * NY; i++)
  {
    for(int j = 0; j <= 2 * NX; j++)
    {
      BioOutPut(i,j) = Bio_Data[j + i * (2 * NX + 1)];
    }
  }

}