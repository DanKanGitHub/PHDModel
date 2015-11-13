#include "readInVelPre.h"

void readInVelPre(int NX,
		  double *data, 
		  char *hVelFilename,
		  char *vVelFilename,
		  char *preFilename) {
  
  int numNodes = (2 * NX + 1) * (2 * NX + 1);

  std::ifstream hVelFile;
  std::ifstream vVelFile;
  std::ifstream preFile;
  
  hVelFile.open(hVelFilename, ios::in);
  vVelFile.open(vVelFilename, ios::in);
  preFile.open(preFilename, ios::in);

  for(int i = 0; i < numNodes; i++) {
    hVelFile >> data[i];
  }
  
  for(int i = 0; i < numNodes; i++) {
    vVelFile >> data[i + numNodes];
  }
  
  for(int i = 0; i < (NX + 1) * (NX + 1); i++) {
    preFile >> data[i + 2 * numNodes];
  }
  
  hVelFile.close();
  vVelFile.close();
  preFile.close();
  
}