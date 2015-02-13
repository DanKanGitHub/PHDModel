#ifndef WriteBioData_H_Guard
#define WriteBioData_H_Guard

#include <fstream>
#include <ostream>

using namespace std;

void WriteBioData(char* filename, 
		  double* Data, 
		  int Size);

#endif