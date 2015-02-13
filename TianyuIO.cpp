#include "DataStruct.h"
#include "MatVec.h"

#define PI (3.14159265359)

//************************************************************************************************************************

//function to output the non-dimensionalized parameters to the screen
void Show_NonDimension_Para() {

  cout << endl << endl <<  " ************* The non-dimensionalized global parameters are ************* " << endl;
  cout << endl << " Vol_Penal = " << Vol_Penal
       << endl << " Gamma_0 = " << Gamma_0
       << endl << " Ca_phase = " << Ca_phase
       << endl << " Gamma_1 = " << Gamma_1
       << endl << " Gamma_2 = " << Gamma_2
       << endl << " Kai_CH  = " << Kai_CH
       << endl << " Lambda_c  = " << Lambda_c
       << endl << " Lambda_b  = " << Lambda_b      
       << endl << " Lambda_s  = " << Lambda_s
       << endl << " Epsilon = " << sub_epsilon
       << endl << " sub_mu  = " << sub_mu
       << endl << " sub_K_c = " << sub_K_c
       << endl << " sub_A   = " << sub_A
       << endl << " a_EPU = " << a_EPU
       << endl << " k_ur   = " << k_ur
       << endl << " k_p   = " <<  k_p 
       << endl << " S_crit = " << S_crit
       << endl << " k_acid_1 = " << K_acid_1
       << endl << " k_acid_2 = " << K_acid_2
       << endl << " k_nh = " << K_nh     
       << endl << " k_w = " << K_w
       << endl << " k_so = " << K_so
       << endl << " Ds  = " << Ds
       << endl << " D_ur  = " << D_ur
       << endl << " D_nh3  = " << D_nh3
       << endl << " D_nh4  = " << D_nh4
       << endl << " D_h2co3  = " << D_h2co3
       << endl << " D_hco3  = " << D_hco3
       << endl << " D_co3  = " << D_co3
       << endl << " D_ca      = " << D_ca
       << endl << " D_cl      = " << D_cl
       << endl << " eta_cac = " << eta_cac      
       << endl << " eta_bio = " << eta_bio
       << endl << " eta_sol = " << eta_sol
       << endl << " rho_cac = " << rho_cac
       << endl << " rho_bio = " << rho_bio
       << endl << " rho_sol = " << rho_sol;


  cout << endl << endl;
}

//************************************************************************************************************************
//function to save the non-dimensionalized parameters to a file
void Save_NonDimension_Para(char* filename){

  ofstream ofile(filename, ios::out);
 
  ofile.flags(ios::scientific | ios::showpos | ios::showpoint );
  ofile.precision(16);   

  if(!ofile)
    {
        cout << endl << "can not open output file " << filename << ", ERROR in Save_NonDimension_Para() at line " << __LINE__ << "of file " << __FILE__ << endl;
      exit(-1);
    }

  ofile << endl << endl <<  " ************* The non-dimensionalized global parameters are ************* " << endl;

  ofile << endl << " Vol_Penal = " << Vol_Penal
       << endl << " Gamma_0 = " << Gamma_0
       << endl << " Ca_phase = " << Ca_phase
       << endl << " Gamma_1 = " << Gamma_1
       << endl << " Gamma_2 = " << Gamma_2
       << endl << " Kai_CH  = " << Kai_CH
       << endl << " Lambda_c  = " << Lambda_c
       << endl << " Lambda_b  = " << Lambda_b      
       << endl << " Lambda_s  = " << Lambda_s
       << endl << " Epsilon = " << sub_epsilon
       << endl << " sub_mu  = " << sub_mu
       << endl << " sub_K_c = " << sub_K_c
       << endl << " sub_A   = " << sub_A
       << endl << " a_EPU = " << a_EPU
       << endl << " k_ur   = " << k_ur
       << endl << " k_p   = " <<  k_p  
       << endl << " S_crit = " << S_crit
       << endl << " k_acid_1 = " << K_acid_1
       << endl << " k_acid_2 = " << K_acid_2
       << endl << " k_nh = " << K_nh     
       << endl << " k_w = " << K_w
       << endl << " k_so = " << K_so    
       << endl << " Ds      = " << Ds
       << endl << " D_ur      = " << D_ur
       << endl << " D_nh3  = " << D_nh3
       << endl << " D_nh4  = " << D_nh4
       << endl << " D_h2co3  = " << D_h2co3
       << endl << " D_hco3  = " << D_hco3
       << endl << " D_co3  = " << D_co3
       << endl << " D_ca      = " << D_ca
       << endl << " D_cl      = " << D_cl
       << endl << " eta_cac = " << eta_cac      
       << endl << " eta_bio = " << eta_bio
       << endl << " eta_sol = " << eta_sol
       << endl << " rho_cac = " << rho_cac
       << endl << " rho_bio = " << rho_bio
       << endl << " rho_sol = " << rho_sol;


  ofile.close();

}    
//************************************************************************************************************************

//function to write a vector into a file as a matrix, used for visualization purpose
//i_max is the number of columns
void Write_Vector(char* filename, const valarray<double>& v, int i_max){

  ofstream ofile(filename, ios::out);
 
  ofile.flags(ios::scientific | ios::showpos | ios::showpoint );
  ofile.precision(16);   

  if(!ofile)
    {
      cout << endl << "can not open output file " << filename << ", ERROR in Write_Vector() in file " << __FILE__ <<" !!! " << endl;
      exit(-1);
    }

  int i, j, j_max;
  j_max = v.size()/i_max;

  for(j = 0; j < j_max; j++)
    {
      for(i = 0; i < i_max-1; i++)
	ofile << v[j*i_max + i] << "  ";
      ofile << v[j*i_max + i_max-1];
      ofile << endl;
    }

  ofile.close();
}

//************************************************************************************************************************

void Write_Array(char* filename, double* v, int length, int i_max){

  ofstream ofile(filename, ios::out);
 
  ofile.flags(ios::scientific | ios::showpos | ios::showpoint );
  ofile.precision(16);   

  if(!ofile)
    {
      cout << endl << "can not open output file " << filename << ", ERROR in Write_Array() in file " << __FILE__ << " !!! " << endl;
      exit(-1);
    }

  int i, j, j_max;
  j_max = length/i_max;

  for(j = 0; j < j_max; j++)
    {
      for(i = 0; i < i_max-1; i++)
	ofile << v[j*i_max + i] << "  ";
      ofile << v[j*i_max + i_max-1];
      ofile << endl;
    }

  ofile.close();
}

//************************************************************************************************************************

//function to read a vector from a file stored as a matrix, used for visualization purpose
//i_max is the number of columns
void Read_Vector(char* filename, valarray<double>& v, int i_max){

  ifstream ifile(filename, ios::in);
 
  if(!ifile)
    {
      cout << endl << "can not open input file " << filename << ", ERROR in Read_Vector() in file " << __FILE__ << " !!! " << endl;
      exit(-1);
    }

  int i, j, j_max;
  j_max = v.size()/i_max;

  for(j = 0; j < j_max; j++)
    {
      for(i = 0; i < i_max; i++)
	ifile >> v[j*i_max + i];
    }

  ifile.close();
}
//************************************************************************************************************************

//function to write a Node_Data variale into a file
void Write_Node_Data(char* filename, const Node_Data& u){

  ofstream ofile(filename, ios::out);
 
  ofile.flags(ios::scientific | ios::showpos | ios::showpoint );
  ofile.precision(16);   
 
  if(!ofile)
    {
      cout << endl << "can not open output file " << filename << ", ERROR in Write_Node_Data() in file" << __FILE__ << "  !!! " << endl;
      exit(-1);
    }

  int i, j, Nx, Ny;
  Nx = u.get_Nx();
  Ny = u.get_Ny();

  for(j = 0; j <= Ny; j++)
    {
      for(i = 0; i < Nx + 1; i++)
	ofile << u[i][j] << "  ";
      ofile << u[Nx+1][j];
      ofile << endl;
    }

  ofile << endl;
  ofile.close();
}

//************************************************************************************************************************
void Write_Node_Data_x_Ave_Ext_N(char* filename, const Node_Data_x_Ave_Ext_N& u){


  ofstream ofile(filename, ios::out);
 
  ofile.flags(ios::scientific | ios::showpos | ios::showpoint );
  ofile.precision(16);   
 
  if(!ofile)
    {
      cout << endl << "can not open output file " << filename << ", ERROR in Write_Node_Data() in file " << __FILE__ << " !!! " << endl;
      exit(-1);
    }

  int i, j, Nx, Ny;
  Nx = u.get_Nx();
  Ny = u.get_Ny();

  for(j = 0; j <= Ny; j++)
    {
      for(i = 0; i < Nx + 2; i++)
	ofile << u[i][j] << "  ";
      ofile << u[Nx+2][j];
      ofile << endl;
    }

  ofile << endl;
  ofile.close();
}


//************************************************************************************************************************
void Write_Node_Data_y_Ave_Ext_N(char* filename, const Node_Data_y_Ave_Ext_N& u){


  ofstream ofile(filename, ios::out);
 
  ofile.flags(ios::scientific | ios::showpos | ios::showpoint );
  ofile.precision(16);   
 
  if(!ofile)
    {
      cout << endl << "can not open output file " << filename << ", ERROR in Write_Node_Data() in file " << __FILE__ << " !!! " << endl;
      exit(-1);
    }

  int i, j, Nx, Ny;
  Nx = u.get_Nx();
  Ny = u.get_Ny();

  for(j = 0; j <= Ny+1; j++)
    {
      for(i = 0; i < Nx + 1; i++)
	ofile << u[i][j] << "  ";
      ofile << u[Nx+1][j];
      ofile << endl;
    }

  ofile << endl;
  ofile.close();
}
//************************************************************************************************************************

//function to write a Node_Data variale into a file
void Write_MultipleRHS(char* filename, const MultipleRHS& u){

  ofstream ofile(filename, ios::out);
 
  ofile.flags(ios::scientific | ios::showpos | ios::showpoint );
  ofile.precision(16);   
 
  if(!ofile)
    {
      cout << endl << "can not open output file " << filename << ", ERROR in Write_MultipleRHS() in file " << __FILE__ << " !!!" << endl;
      exit(-1);
    }

  int i, j, dim, col_num;
  col_num = u.get_col_num();
  dim = u.get_dim();

  for(j = 0; j < col_num; j++)
    {
      for(i = 0; i < dim-1; i++)
	ofile << u[j][i] << "  ";
      ofile << u[j][dim-1];
      ofile << endl;
    }

  ofile << endl;
  ofile.close();
}


//************************************************************************************************************************
void Write_All_Node_Data(char* SaveDir, int file_NO, const Node_Data& u, const Node_Data& v, const Node_Data& phi_c, const Node_Data& phi_b, const Node_Data& Ur, 
                         const Node_Data& NH3, const Node_Data& NH4, const Node_Data& H2CO3, const Node_Data& HCO3, const Node_Data& CO3, const Node_Data& Ca,
                         const Node_Data& Cl, const Node_Data& H, const Node_Data& OH, const Node_Data& NT, const Node_Data& CT, const Node_Data& EP, 
                         const Node_Data& SS, const Node_Data& c, double T_current,  const Control_Parameter& CP){

  char debfile[80];

  sprintf(debfile,"%su_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,u);
 
  sprintf(debfile,"%sv_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,v);

  sprintf(debfile,"%sphic_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,phi_c);

  sprintf(debfile,"%sphib_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,phi_b);      

  sprintf(debfile,"%sUr_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,Ur);

  sprintf(debfile,"%sNH3_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,NH3);

  sprintf(debfile,"%sNH4_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,NH4);

  sprintf(debfile,"%sH2CO3_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,H2CO3);

  sprintf(debfile,"%sHCO3_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,HCO3);

  sprintf(debfile,"%sCO3_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,CO3);

  sprintf(debfile,"%sCa_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,Ca);

  sprintf(debfile,"%sCl_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,Cl);

  sprintf(debfile,"%sH_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,H);  

  sprintf(debfile,"%sOH_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,OH);  

//   sprintf(debfile,"%sNT_%d",SaveDir,file_NO);
//   Write_Node_Data(debfile,NT);  

//   sprintf(debfile,"%sCT_%d",SaveDir,file_NO);
//   Write_Node_Data(debfile,CT);
  

  if(CP.EP_flag == 1)
    {
      sprintf(debfile,"%sEP_%d",SaveDir,file_NO);
      Write_Node_Data(debfile,EP);  
    }

  sprintf(debfile,"%sS_%d",SaveDir,file_NO);
  Write_Node_Data(debfile,SS);

//   if(T_current < CP.cutoff_c_time)
//     {
//       sprintf(debfile,"%sc_%d",SaveDir,file_NO);
//       Write_Node_Data(debfile,c);
//     }
}


//************************************************************************************************************************
void Screen_Display(int n, double tdiff, double dt_n, double ave_viscosity, double t_curr, double dx, double dy, const Node_Data& u, const Node_Data& v, 
                    const Node_Data& phi_c, const Node_Data& phi_b, const Node_Data& Ur, const Node_Data& NH3, const Node_Data& NH4, const Node_Data& H2CO3, 
                    const Node_Data& HCO3, const Node_Data& CO3, const Node_Data& Ca, const Node_Data& Cl, const Node_Data& H, const Node_Data& OH, 
                    const Node_Data& SS, const Node_Data& c){

  cout << endl  << "_______________________________________________________________________________________" 
       << endl  <<  endl << " At time t = " << t_curr << " after " << n << " time steps "
       << endl  << " The clock time is: ";
 
  system("date");

  cout << endl << endl << " Time step " << n << " took " << tdiff << " seconds, dt_n = " << dt_n 
       << ", eta_ave = " <<  ave_viscosity << endl 
       << " L_infty norm of u is " << u.L_inf_norm() << ", L_infty norm of v is "  << v.L_inf_norm()
       << endl << " max(phi_c) = " << phi_c.L_inf_norm() << ", total(phi_c) = " << phi_c.sum()*dx*dy
       << endl << " max(phi_b) = " << phi_b.L_inf_norm() << ", total(phi_b) = " << phi_b.sum()*dx*dy
       << endl << " max(Ur) = " << Ur.L_inf_norm() << ", total(Ur) = " << Ur.sum()*dx*dy
       << endl << " max(NH3) = " << NH3.L_inf_norm() << ", total(NH3) = " << NH3.sum()*dx*dy
       << endl << " max(NH4) = " << NH4.L_inf_norm() << ", total(NH4) = " << NH4.sum()*dx*dy
       << endl << " max(H2CO3) = " << H2CO3.L_inf_norm() << ", total(H2CO3) = " << H2CO3.sum()*dx*dy
       << endl << " max(HCO3) = " << HCO3.L_inf_norm() << ", total(HCO3) = " << HCO3.sum()*dx*dy
       << endl << " max(CO3) = " << CO3.L_inf_norm() << ", total(CO3) = " << CO3.sum()*dx*dy
       << endl << " max(Ca) = " << Ca.L_inf_norm() << ", total(Ca) = " << Ca.sum()*dx*dy
       << endl << " max(Cl) = " << Cl.L_inf_norm() << ", total(Cl) = " << Cl.sum()*dx*dy
       << endl << " max(H) = " << H.L_inf_norm() << ", total(H) = " << H.sum()*dx*dy
       << endl << " max(OH) = " << OH.L_inf_norm() << ", total(OH) = " << OH.sum()*dx*dy
       << endl << " max(S) = " << SS.L_inf_norm() << ", total(S) = " << SS.sum()*dx*dy
       << endl << " max(c) = " << c.L_inf_norm() << ", total(c) = " << c.sum()*dx*dy << endl;
}

//************************************************************************************************************************

//overload the output operator << for valarray<double>
ostream& operator<<(ostream& _ostr, const valarray<double>& u){

  int i;

  std::ios::fmtflags oldflags = _ostr.flags(ios::scientific | ios::showpos | ios::showpoint );

  int oldprecision =  _ostr.precision(16); 

  for(i = 0; i < u.size(); i++)
    _ostr << u[i] << "  ";

  _ostr << endl << endl;

  _ostr.precision(oldprecision);

  _ostr.flags(oldflags);

  return _ostr;
}

//overload the input operator << for valarray<double>
istream& operator>>(istream& _istr, valarray<double>& u){

  int i;

  for(i = 0; i < u.size(); i++)
    _istr >> u[i];

  return _istr;
}
//************************************************************************************************************************
