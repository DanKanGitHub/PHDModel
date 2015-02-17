// Writes data to filename using Tianyu's program for inspiration
// Size = NX

#include "WriteVelPreData.h"

void WriteVelPreData(char* HorVelfilename, 
		     char* VertVelfilename,
		     char* Prefilename, 
		     double* Data, 
		     int Size)
{
  
  std::ofstream HorVelofile(HorVelfilename, std::ios::out);
  
  HorVelofile.flags(std::ios::scientific | std::ios::showpos | std::ios::showpoint );
  HorVelofile.precision(16);
  
  ofstream VertVelofile(VertVelfilename, std::ios::out);
  
  VertVelofile.flags(std::ios::scientific | std::ios::showpos | std::ios::showpoint );
  VertVelofile.precision(16);
  
  for(int i = 0; i < 2 * Size + 1; i++)
  {
    
    for(int j = 0; j < 2 * Size + 1; j++)
    {
      
      HorVelofile << Data[i * (2 * Size + 1) + j] << "  ";
      VertVelofile << Data[i * (2 * Size + 1) + j + (2 * Size + 1) * (2 * Size + 1)] << "  ";
      
    }

    HorVelofile << endl;

    VertVelofile << endl;
    
  }
  
  std::ofstream Preofile(Prefilename, std::ios::out);
  
  Preofile.flags(std::ios::scientific | std::ios::showpos | std::ios::showpoint );
  Preofile.precision(16);
  
  for(int i = 0; i < Size + 1; i++)
  {
    
    for(int j = 0; j < Size + 1; j++)
    {
      
      //(2 * Size + 1) * (2 * Size + 1) is the number of Hor Vel terms and the number of Vert Vel terms
      Preofile << Data[i * (Size + 1) + j + 2 * (2 * Size + 1) * (2 * Size + 1)] << "  ";
      
    }

    Preofile << endl;

  }
  
  HorVelofile.close();
  VertVelofile.close();
  Preofile.close();
  
}