#ifndef QuadPtWt_H_Guard
#define QuadPtWt_H_Guard

#include <cmath>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_SerialDenseVector E_SDV;

void QuadPtWt(int NTRIQUAD, 
	      E_SDM & Tri_Quad_Pt, 
	      E_SDV & Tri_Quad_Wt);

#endif