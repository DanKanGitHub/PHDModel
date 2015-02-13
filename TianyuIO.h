#ifndef _IO_H
#define _IO_H

#include "Build.h"

using namespace std;

//function to output the non-dimensionalized parameters to the screen
void Show_NonDimension_Para();

//function to save the non-dimensionalized parameters to a file
void Save_NonDimension_Para(char* filename);

//function to write a vector into a file as a matrix, used for visualization purpose
//i_max is the number of columns
void Write_Vector(char* filename, const valarray<double>& v, int i_max);

//function to write an array into a file as a matrix, used for visualization purpose
//i_max is the number of columns
void Write_Array(char* filename, double* v, int length, int i_max);

//************************************************************************************************************************

//function to read a vector from a file stored as a matrix, used for visualization purpose
//i_max is the number of columns
void Read_Vector(char* filename, valarray<double>& v, int i_max);

//function to write a Node_Data variale into a file
void Write_Node_Data(char* filename, const Node_Data& u);

//function to write a Node_Data_x_Ave_Ext_N variale into a file
void Write_Node_Data_x_Ave_Ext_N(char* filename, const Node_Data_x_Ave_Ext_N& u);

//function to write a Node_Data_x_Ave_Ext_N variale into a file
void Write_Node_Data_y_Ave_Ext_N(char* filename, const Node_Data_y_Ave_Ext_N& u);

//function to write a Node_Data variale into a file
void Write_MultipleRHS(char* filename, const MultipleRHS& u);

//function to write all important Node_Data variables to file
void Write_All_Node_Data(char* SaveDir, int file_NO, const Node_Data& u, const Node_Data& v, const Node_Data& phi_c, const Node_Data& phi_b, const Node_Data& Ur, 
                         const Node_Data& NH3, const Node_Data& NH4, const Node_Data& H2CO3, const Node_Data& HCO3, const Node_Data& CO3, const Node_Data& Ca,
                         const Node_Data& Cl, const Node_Data& H, const Node_Data& OH, const Node_Data& NT, const Node_Data& CT, const Node_Data& EP, 
                         const Node_Data& SS, const Node_Data& c, double T_current,  const Control_Parameter& CP);

//function to display information on screen 
void Screen_Display(int n, double tdiff, double dt_n, double ave_viscosity, double t_curr, double dx, double dy, const Node_Data& u, const Node_Data& v, 
                    const Node_Data& phi_c, const Node_Data& phi_b, const Node_Data& Ur, const Node_Data& NH3, const Node_Data& NH4, const Node_Data& H2CO3, 
                    const Node_Data& HCO3, const Node_Data& CO3, const Node_Data& Ca, const Node_Data& Cl, const Node_Data& H, const Node_Data& OH, 
                    const Node_Data& SS, const Node_Data& c);


//************************************************************************************************************************

//overload the output operator << for valarray<double>
ostream& operator<<(ostream& _ostr, const valarray<double>& u);

//overload the input operator << for valarray<double>
istream& operator>>(istream& _istr, valarray<double>& u);



#endif
