#ifndef WriteVelPreData_H_Guard
#define WriteVelPreData_H_Guard

#include <fstream>
#include <iostream>

using namespace std;

void WriteVelPreData(char* HorVelfilename, 
		     char* VertVelfilename,
		     char* Prefilename, 
		     double* Data, 
		     int Size);

#endif