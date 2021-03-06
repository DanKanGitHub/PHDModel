// Writes data to filename using Tianyu's program for inspiration
// Size = NX + 1 for Stress
// Size = 2 * NX + 1 for Bio-Film

#include "WriteStrData.h"

void WriteStrData(char* XXfilename,
		  char* XYfilename,
		  char* YYfilename,
		  E_SDM Data, 
		  int Size)
{
  
  std::ofstream XXofile(XXfilename, ios::out);
  
  XXofile.flags(ios::scientific | ios::showpos | ios::showpoint );
  XXofile.precision(16);
  
  std::ofstream XYofile(XYfilename, ios::out);
  
  XYofile.flags(ios::scientific | ios::showpos | ios::showpoint );
  XYofile.precision(16);
  
  std::ofstream YYofile(YYfilename, ios::out);
  
  YYofile.flags(ios::scientific | ios::showpos | ios::showpoint );
  YYofile.precision(16);

  
  for(int i = 0; i < Size; i++)
  {
    
    for(int j = 0; j < Size; j++)
    {
      
      XXofile << Data(i * Size + j,0) << "  ";
      XYofile << Data(i * Size + j,1) << "  ";
      YYofile << Data(i * Size + j,2) << "  ";
      
    }
    
//     XXofile << Data(i * Size + Size - 1,0);
    XXofile << endl;
    
//     XYofile << Data(i * Size + Size - 1,1);
    XYofile << endl;
    
//     YYofile << Data(i * Size + Size - 1,2);
    YYofile << endl;
    
  }
  
  XXofile.close();
  XYofile.close();
  YYofile.close();
  
}