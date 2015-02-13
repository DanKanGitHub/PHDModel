#ifndef WriteStrData_H_Guard
#define WriteStrData_H_Guard

#include <fstream>
#include <ostream>

#include "Epetra_SerialDenseMatrix.h"

typedef Epetra_SerialDenseMatrix E_SDM;

using namespace std;

void WriteStrData(char* XXfilename,
		  char* XYfilename,
		  char* YYfilename,
		  E_SDM Data, 
		  int Size);

#endif