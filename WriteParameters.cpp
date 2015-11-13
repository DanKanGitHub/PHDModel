#include "WriteParameters.h"

void WriteParameters(char* filename,
		     int NX,
		     int NY,
		     int NGP,
		     int N_TRI_QUAD,
		     int MAX_TIME_STEP_NUM,
		     int Write_Time_Steps_Skipped,
		     double SOL_NEW_VIS,
		     double POL_NEW_VIS,
		     double RELAX_TIME,
		     double RETARD_TIME,
		     double POL_DENSITY,
		     double SOL_DENSITY,
		     double T_SCALE,
		     double L_SCALE,
		     double DIFF_COEFF,
		     double BIO_DIFF_COEFF,
		     double TIMESTEPSCALAR)
{
  
  std::ofstream ofile(filename, ios::out);
  
  ofile.flags(ios::scientific | ios::showpos | ios::showpoint );
  ofile.precision(16);

  std::ofstream ParameterListOutFileName(filename, std::ios::out);
  
  ParameterListOutFileName << "Number of slices in the horizontal direction." << endl;
  ParameterListOutFileName << "NX = " << NX << endl;
  
  ParameterListOutFileName << "Number of slices in the vertical direction." << endl;
  ParameterListOutFileName << "NY = " << NY << endl;
  
  ParameterListOutFileName << "Number of Gauss points used when integrating on the boundary for natural boundary conditions" << endl;
  ParameterListOutFileName << "NGP = " << NGP << endl;
  
  ParameterListOutFileName << "Number of points used in Gaussian quadrature inside of elements." << endl;
  ParameterListOutFileName << "N_TRI_QUAD = " << N_TRI_QUAD << endl;
  
  ParameterListOutFileName << "Maximum number of time steps that can be taken." << endl;
  ParameterListOutFileName << "MAX_TIME_STEP_NUM = " << MAX_TIME_STEP_NUM << endl;
  
  ParameterListOutFileName << "Number of time steps that are skipped between writing data to files." << endl;
  ParameterListOutFileName << "Write_Time_Steps_Skipped = " << Write_Time_Steps_Skipped << endl;

  ParameterListOutFileName << "Solvent Newtonian Viscosity." << endl;
  ParameterListOutFileName << "SOL_NEW_VIS = " << SOL_NEW_VIS << endl;
  
  ParameterListOutFileName << "Polymer Newtonian Viscosity." << endl;
  ParameterListOutFileName << "POL_NEW_VIS = " << POL_NEW_VIS << endl;

  ParameterListOutFileName << "Relaxation Time." << endl;
  ParameterListOutFileName << "RELAX_TIME = " << RELAX_TIME << endl;
  
  ParameterListOutFileName << "Retardation Time." << endl;
  ParameterListOutFileName << "RETARD_TIME = " << RETARD_TIME << endl;
  
  ParameterListOutFileName << "Density of the polymer." << endl;
  ParameterListOutFileName << "POL_DENSITY = " << POL_DENSITY << endl;
  
  ParameterListOutFileName << "Density of the solvent." << endl;
  ParameterListOutFileName << "SOL_DENSITY = " << SOL_DENSITY << endl;
  
  ParameterListOutFileName << "The time scale of the model, seconds." << endl;
  ParameterListOutFileName << "T_SCALE = " << T_SCALE << endl;
  
  ParameterListOutFileName << "The length scale of the model, meters." << endl;
  ParameterListOutFileName << "L_SCALE = " << L_SCALE << endl;
  
  ParameterListOutFileName << "Diffusion Coefficient for advection diffusion model." << endl;
  ParameterListOutFileName << "DIFF_COEFF = " << DIFF_COEFF << endl;
  
  ParameterListOutFileName << "Diffusion Coefficient for the softening of polymer initial condition." << endl;
  ParameterListOutFileName << "BIO_DIFF_COEFF = " << BIO_DIFF_COEFF << endl;
  
  ParameterListOutFileName << "The scalar which converts the spatial step size to the time step size." << endl;
  ParameterListOutFileName << "TIMESTEPSCALAR = " << TIMESTEPSCALAR << endl;
  
  ofile.close();
  
}