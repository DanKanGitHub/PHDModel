#include "readInStr.h"

void readInStr(int NX,
		  E_SDM data, 
		  char *XXFilename,
		  char *XYFilename,
		  char *YYFilename) {
  
  int numNodes = (2 * NX + 1) * (2 * NX + 1);

  std::ifstream XXFile;
  std::ifstream XYFile;
  std::ifstream YYFile;
  
  XXFile.open(XXFilename, ios::in);
  XYFile.open(XYFilename, ios::in);
  YYFile.open(YYFilename, ios::in);

  for(int i = 0; i < numNodes; i++) {
    XXFile >> data(i,0);
    XYFile >> data(i,1);
    YYFile >> data(i,2);
  }

  XXFile.close();
  XYFile.close();
  YYFile.close();
  
}