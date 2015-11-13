#include "readInBio.h"

void readInBio(int NX,
	       double *Bio, 
	       char *fileName) {
  
  int numNodes = (2 * NX + 1) * (2 * NX + 1);

  std::ifstream iBioFile;
  
  iBioFile.open(fileName, ios::in);

  for(int i = 0; i < numNodes; i++) {
    iBioFile >> Bio[i];
  }
  
  iBioFile.close();
  
}