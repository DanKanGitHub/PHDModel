// Viscosity * 0.5 * (Gradu + GraduTrans)
// The other values not discribed here are initialized to zero.
// While Str is symmetric all code takes advatage of this property so it isn't necesary to actually code
// more than one term since the main diagonal of the stress tensor is initiially zero.

#include "InitialStr.h"

void InitialStr(int Nem, 
		int Pre_Npe,
		double Sol_Vis,
		double Poly_Vis,
		double Sol_Density,
		double Poly_Density,
		double L_ZERO,
		double T_ZERO,
		E_ISDM Vel_Nod,
		E_ISDM Pre_Nod,
		E_SDM Glxy,
		double *Bio,
		E_SDM & Str) 
{
  double y, Effective_Vis, Bio_Weight, RetardDivRelax, U_Zero, P_Zero, Effective_Den, Effective_Pol_Vis;
  int Bio_GP, Str_GP, Count;
  
//   for (int i = 0; i <= Nnm-1; i++)
//   {
//     y = Glxy(i,1);
//     Str(i, 1) = NewVis * (3.0 - 6.0 * y);
//   }

  U_Zero = L_ZERO / T_ZERO;

  for(int Ne = 0; Ne < Nem; Ne++)
  {
    for(int Np = 0; Np < Pre_Npe; Np++)
    {

      Bio_GP = Vel_Nod(Ne, Np) - 1;
      Str_GP = Pre_Nod(Ne, Np) - 1; // Global Point
      
//       std::cout << "Str(Gp,1) = " << Str(Gp,1) << endl;
      
      if(Str(Str_GP,1) == 1000000)
      {
	
// 	std::cout << "Made it!" << endl;
	
	// Returns the bio-film volume fraction.
	Bio_Weight = BioWeightFunc(Bio[Bio_GP]);

	RetardDivRelax = RetardationDividedByRelaxation(Bio_Weight);

	Effective_Den = (1 - Bio_Weight) * Sol_Density + Bio_Weight * Poly_Density;
	
	P_Zero = Effective_Den * U_Zero * U_Zero;
	
	Effective_Vis = (Sol_Vis * (1 - Bio_Weight) + Poly_Vis * Bio_Weight) / (T_ZERO * P_Zero);
	
	Effective_Pol_Vis = (1 - RetardDivRelax) * Effective_Vis;

	y = Glxy(Str_GP,1);
	
	Str(Str_GP, 1) = Effective_Pol_Vis * (6.0 - 12.0 * y);
	
      }
      
    }
  }

}

// GaussPt_Bio_Weight = BioWeightFunc(GaussPt_Bio);
// Effective_Vis = Sol_Vis + Poly_Vis * GaussPt_Bio_Weight;