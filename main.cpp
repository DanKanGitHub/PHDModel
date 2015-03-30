// Main program for my 2D FEM model of viscco-elastic fluid flow.
// Used Prathish's code for MPI implementation and IfPack

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "AztecOO_Operator.h"
#include "Epetra_Import.h"

#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"

#include <fstream>
#include <ostream>
#include "VelBioMesh2d.h"
#include "PreMesh2d.h"
#include "StrMesh2d.h"
#include "SparseAssembly.h"
#include "EquibSparseAssembly.h"
#include "QuadPtWt.h"
#include "InitialVel.h"
#include "InitialStr.h"
#include "InitialBio.h"
#include "ElementNeigh.h"
#include "Conformation.h"
#include "InitMatZero.h"
#include "ProcNodePartitionClass.h"
#include "InitMatZeroBio.h"
#include "BioSparseAssembly.h"
#include "BioDiffSparseAssembly.h"
#include "WriteVelPreData.h"
#include "WriteBioData.h"
#include "WriteStrData.h"
#include "VectorDiffNorm.h"
#include "VectorNorm.h"
#include "InitialGuess.h"
#include "GaussDepartureFeet.h"
#include "StrNodeDepartureFeet.h"
#include "StressReassemble.h"

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_BlockRelaxation.h"
#include "Ifpack_SparseContainer.h"
#include "Ifpack_Amesos.h"
#include "Ifpack_Graph.h"
#include "Ifpack_Graph_Epetra_CrsGraph.h"
#include "Ifpack_Graph_Epetra_RowMatrix.h"
#include "Ifpack_DenseContainer.h"
#include "Ifpack.h"
#include <Ifpack_Preconditioner.h>

typedef Epetra_CrsMatrix E_CM;
typedef Epetra_Map E_Mp;
typedef Epetra_Vector E_V;

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_IntSerialDenseMatrix E_ISDM;
typedef Epetra_SerialDenseVector E_SDV;
typedef Epetra_IntSerialDenseVector E_ISDV;

using namespace std;

//Constants
const int NX = 32;			// Number of element intervals in the horizontal direction
const int NY = 32;
const int NGP = 4;			// Number of Gauss points in numerical quadrature, used on the boundary
const int N_TRI_QUAD = 7;		// Number of Gauss points in numerical quadrature, used in the element
const int MAX_TIME_STEP_NUM = 1000;	// Maximum number of time interations

const int VEL_FLAG = 2;			// Program flags 1 = linear, 2 = quadratic
const int PRE_FLAG = 1;
const int STRESS_FLAG = 2;

const int Write_Time_Steps_Skipped = 1;

const double XL = 0.0;			// coordinate of left boundary element
const double XR = 1.0;			// coordinate of right boundary element
const double YB = 0.0;			// Coordinate of the bottom boundary of the domain.
const double YT = 1.0;			// Coordinate of the top boundary of the domain.
// Tianyu's Phasse field paper
const double SOL_NEW_VIS = 2.0;	// Solvent Newtonian Viscosity, Water dynamic viscosity at 25 C.
const double POL_NEW_VIS = 100.0; //10000.0;	// Polymer Newtonian Viscosity, Honey dynamic viscosity at 25 C.

// Issac's Viscoelastic Fluid Description paper
const double RELAX_TIME = 1080.0; //100.0;	// Relaxation time.
const double RETARD_TIME = 1079.0; //10.0;	// Retardation time, !!!must be less than Relaxation time!!!

const double DENSITY_BIO = 1.0;	// Density Biofilm, kg/m^3 at 25 C.
const double DENSITY_WATER = 1.0;	// Density Water, kg/m^3 at 25 C.
const double T_ZERO = 60.0;		// seconds (I read Issac's paper for conformation.)
const double L_ZERO = 0.01;		// meters
const double TOL = 0.000001;
const double BIO_TOL = 0.0000001;
const double DIFF_TOL = 0.000000000001;		// 10^(-12)
// Hits 0.0000001 for 64X64!
const double EQUIB_TOL = 0.00001;		// The tolerance that determines if the pre-growth simulation has reached equilibrium
const double DIFF_COEFF = 0.01; // 0.00001;	// This is a guess.  Diffusion coefffcient for eps diffusion into water
const double BIO_DIFF_COEFF = 1.0;
const double TIMESTEPSCALAR = 0.01;

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  int NumProc;
  int myid;
  MPI_Comm_size (MPI_COMM_WORLD, &NumProc);
  MPI_Comm_rank (MPI_COMM_WORLD, &myid);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  
  // variable decalrations
  // Velocity Related
  // Velocity is NNMX2, Depart_Vel is Gaus_Pt_NumX2
  E_SDM Vel_Glxy, Vel_Old;
  E_ISDM Vel_Nod, Vel_Nod_BC_Hor, Vel_Nod_BC_Ver;
  E_SDV Vel_Norm;
  int Vel_Npe, Vel_Nnm;
  double Tol;
  
  // Presure Related
  E_SDM Pre_Glxy;
  E_ISDM Pre_Nod, Pre_Nod_BC;
  E_SDV Pre_Norm;
  int Pre_Npe, Pre_Nnm;

  // Stress Related
  // Stess is NNMX3
  E_SDM Str, Str_New, Str_Old, Str_Glxy;
  E_ISDM Str_Nod;
  int Str_Npe, Str_Write_Size, Str_Nnm;
  
  // Bio Related
  E_SDV Bio_Old, Bio_Norm;
  E_ISDM Bio_Nod_BC;
  
  // Other
  E_SDM GAUSPT, GAUSWT, Tri_Quad_Pt, GaussDepartFootx, GaussDepartFooty, StrNodeDepartFootx, StrNodeDepartFooty;
  E_ISDM Ele_Neigh, GaussDepartElement, StrNodeDepartElement;
  E_SDV Tri_Quad_Wt;
  int Nem, L, n, File_No, Diff_File_No, Bio_File_No;
  double dx, dy, Time_Step, SolnNorm, SolnDiffNorm, EquibTol;
  char HorVelFilename[100];
  char VertVelFilename[100];
  char PreFilename[100];
  char StrXXFilename[100];
  char StrXYFilename[100];
  char StrYYFilename[100];
  char BioFilename[100];

  // Gaussian quadrature points and weights for 1d boundary integration, NGP = 1,2 3, or 4
  // In this code it is intended that NGP = 4
  GAUSPT.Shape(4,4);
  GAUSWT.Shape(4,4);

  GAUSPT(0,1) = -1.0/sqrt(3);
  GAUSPT(0,3) = 0.33998194;
  GAUSPT(1,1) = 1.0/sqrt(3);
  GAUSPT(1,2) = 0.7745966692;
  GAUSPT(1,3) = -0.33998194;
  GAUSPT(2,2) = -0.7745966692;
  GAUSPT(2,3) = 0.86113631;
  GAUSPT(3,3) = -0.86113631;
  
  GAUSWT(0,0) = 2.0;
  GAUSWT(0,1) = 1.0;
  GAUSWT(0,2) = 8.0/9.0;
  GAUSWT(0,3) = 0.65214515;
  GAUSWT(1,1) = 1.0;
  GAUSWT(1,2) = 5.0/9.0;
  GAUSWT(1,3) = 0.65214515;
  GAUSWT(2,2) = 5.0/9.0;
  GAUSWT(2,3) = 0.34785485;
  GAUSWT(3,3) = 0.34785485;

  Nem = 2*NX*NY; // Number of triangular nodes in the mesh.
  
  // Step sizes
  dx = (XR - XL) / NX;	// mesh size in the x direction
  dy = (YT - YB) / NY;	// mesh size in the y direction
  
  // Choose Time_Step so small that even at maximum velocity the  Gauss point closest to the edge of an 
  // element cannot move out of that element in a single time step
  if ( dx <= dy)
  {
    Time_Step = TIMESTEPSCALAR * dx;
  }
  else
  {
    Time_Step = TIMESTEPSCALAR * dy;
  }

  // preprocessor, generate mesh and connectivity matrix
  VelBioMesh2d(dx, 
	      dy,
	      XL,
	      YB,
	      NX, 
	      NY,
	      Nem, 
	      Vel_Nnm, 		// output
	      Vel_Npe, 		// output
	      Vel_Glxy, 	// output
	      Vel_Nod,		// output
	      Vel_Nod_BC_Hor,	// output
	      Vel_Nod_BC_Ver,	// output
	      Bio_Nod_BC);	// output

  PreMesh2d(dx, 
	    dy, 
	    XL,
	    YB,
	    NX, 
	    NY,
	    Nem, 
	    Pre_Nnm, 		// output
	    Pre_Npe, 		// output
	    Pre_Glxy, 		// output
	    Pre_Nod,		// output
	    Pre_Nod_BC);	// output
  
  StrMesh2d(dx, 
	     dy, 
	     XL,
	     YB,
	     STRESS_FLAG,
	     NX, 
	     NY,
	     Nem, 
	     Str_Nnm, 
	     Str_Npe, 
	     Str_Glxy, 
	     Str_Nod);

  GaussDepartFootx.Shape(Nem, N_TRI_QUAD);
  GaussDepartFooty.Shape(Nem, N_TRI_QUAD);
  GaussDepartElement.Shape(Nem, N_TRI_QUAD);
  
  StrNodeDepartFootx.Shape(Nem, Str_Npe);
  StrNodeDepartFooty.Shape(Nem, Str_Npe);
  StrNodeDepartElement.Shape(Nem, Str_Npe);
  
  // Assigns the correct row length for printing to files
  if(STRESS_FLAG == 1)
  {
    Str_Write_Size = NX + 1;
  }
  else
  {
    Str_Write_Size = 2 * NX + 1;
  }
  
  ProcNodePartitionClass Proc_Node_Part;
  
  Proc_Node_Part.ProcNodePartition(myid,
				    NumProc, 
				    Vel_Nnm, 
				    Pre_Nnm, 
				    Vel_Npe,
				    Pre_Npe,
				    Vel_Nod,
				    Pre_Nod,
				    Nem);

  // create a linear map
  int NumGlobalElements_VP = 2 * Vel_Nnm + Pre_Nnm;
  int NumGlobalElements_Bio = Vel_Nnm;

  E_Mp Local_Proc_Map_VP(-1, Proc_Node_Part.Nodes_Per_Proc_VP, Proc_Node_Part.My_Proc_Nodes_VP.Values(), 0, Comm);
  E_Mp Full_Sol_Map_VP(NumGlobalElements_VP, NumGlobalElements_VP, 0, Comm);
  Epetra_Import CompleteSolution_Importer_VP(Full_Sol_Map_VP, Local_Proc_Map_VP);

  E_Mp Local_Proc_Map_Bio(-1, Proc_Node_Part.Nodes_Per_Proc_Bio, Proc_Node_Part.My_Proc_Nodes_Bio.Values(), 0, Comm);
  E_Mp Full_Sol_Map_Bio(NumGlobalElements_Bio, NumGlobalElements_Bio, 0, Comm);
  Epetra_Import CompleteSolution_Importer_Bio(Full_Sol_Map_Bio, Local_Proc_Map_Bio);

  // Create an Epetra_Matrix
  //For my mesh I have, at most, 45 (19 for each component of velocity and 7 from pressure) nonzeros per row
  E_CM A_VP(Copy, Local_Proc_Map_VP, 10 * (2 * Vel_Npe + Pre_Npe));
  E_CM A_Bio(Copy, Local_Proc_Map_Bio, 10 * Vel_Npe);
  
  // Create Solution Vectors
  //Local RHS and soln vectors
  E_V x_VP(Local_Proc_Map_VP), b_VP(Local_Proc_Map_VP);
  E_V x_Bio(Local_Proc_Map_Bio), b_Bio(Local_Proc_Map_Bio);
  
  // Global Soln vectors
  E_V Soln_Cur_t(Full_Sol_Map_VP), Soln_Cur_L(Full_Sol_Map_VP), Soln_Next_L(Full_Sol_Map_VP);
  E_V Bio_Soln_Cur_t(Full_Sol_Map_Bio), Bio_Soln_Next_t(Full_Sol_Map_Bio);

  Soln_Cur_t.PutScalar(0.0);
  
  // Initialize the matrix A to all zeros.  This is Prathish's idea
  InitMatZero(Vel_Npe, 
	      Pre_Npe,
	      Nem, 
	      Vel_Nnm,
	      Pre_Nnm,
	      Vel_Nod, 
	      Pre_Nod,
	      Vel_Nod_BC_Hor,
	      Vel_Nod_BC_Ver,
	      Pre_Nod_BC,
	      Proc_Node_Part.All_Proc_Nodes_VP,
	      Proc_Node_Part.My_Proc_Eles,
	      myid,
	      A_VP);

  InitMatZeroBio(Vel_Npe, 
		  Nem,
		  Vel_Nnm,
		  Vel_Nod, 
		  Bio_Nod_BC,
		  Proc_Node_Part.All_Proc_Nodes_Bio,
		  Proc_Node_Part.My_Proc_Eles,
		  myid,
		  A_Bio);

  A_VP.FillComplete();
  A_Bio.FillComplete();
  
  // Velocity, Pressure and Stress Data
  Str.Shape(Str_Nnm,3); // Only Pre_Nnm X 3 because the stress is symmetric
  Str_New.Shape(Str_Nnm,3); // Only Pre_Nnm X 3 because the stress is symmetric
  Str_New.Shape(Str_Nnm,3); // Only Pre_Nnm X 3 because the stress is symmetric
  Str_Old.Shape(Str_Nnm,3);

  // Initialize the sizes of the LHS matrices and RHS vectors
  Vel_Norm.Size(2 * Vel_Nnm);
  Pre_Norm.Size(Pre_Nnm);
  Bio_Norm.Size(Vel_Nnm);

  // Determine intial surfaces and stress tensor nodal values
  InitialVel(Vel_Nnm, 
	     Vel_Glxy, 
	     Soln_Cur_t.Values()); // Size Vel_NnmX2, output

  // Initialize Biofilm vector
  InitialBio(NX,
	     Vel_Glxy,
	     Bio_Soln_Cur_t.Values());
  
  for(int i = 0; i < Str_Nnm; i++)
  {
    Str(i,1) = 1000000;
  }

  InitialStr(Nem, 
	    Str_Npe,
	    SOL_NEW_VIS,
	    POL_NEW_VIS,
	    DENSITY_WATER,
	    DENSITY_BIO,
	    L_ZERO,
	    T_ZERO,
	    Vel_Nod,
	    Str_Nod,
	    Str_Glxy, 
	    Bio_Soln_Cur_t.Values(),
	    Str); // Size Pre_NnmX3, the stress tensor is symmetric, output

  File_No = 0;
  Diff_File_No = 0;
  Bio_File_No = 0;
  
  sprintf(HorVelFilename, "./Data/Velocity/HorVel/HorVelFile_00%d.data", File_No);
  sprintf(VertVelFilename, "./Data/Velocity/VertVel/VertVelFile_00%d.data", File_No);
  sprintf(PreFilename, "./Data/Pressure/PreFile_00%d.data", File_No);
  
  WriteVelPreData(HorVelFilename, 
		  VertVelFilename,
		  PreFilename , 
		  Soln_Cur_t.Values(), 
		  NX);
  
  sprintf(BioFilename, "./Data/Biofilm/AdvecDiff/BioFile_00%d.data", Bio_File_No);
  
  WriteBioData(BioFilename,
		Bio_Soln_Cur_t.Values(), 
		2 * NX + 1);
  
  sprintf(BioFilename, "./Data/Biofilm/DiffusionOnly/BioFile_00%d.data", Diff_File_No);
  
  WriteBioData(BioFilename,
		Bio_Soln_Cur_t.Values(), 
		2 * NX + 1);
  
  sprintf(StrXXFilename, "./Data/Stress/XX/StrXXFile_00%d.data", File_No);
  sprintf(StrXYFilename, "./Data/Stress/XY/StrXYFile_00%d.data", File_No);
  sprintf(StrYYFilename, "./Data/Stress/YY/StrYYFile_00%d.data", File_No);
  
  WriteStrData(StrXXFilename,
	       StrXYFilename,
	       StrYYFilename,
	       Str, 
	       Str_Write_Size);
  
  File_No++;
  Diff_File_No++;
  Bio_File_No++;

  // There are at most 3 neightbors in my triangular mesh.
  Ele_Neigh.Shape(Nem,3);

  // Determine which global elements share each side of each element
  ElementNeigh(Nem,
	       NX,
	       Ele_Neigh);	// output

  //The points and weights for quadrature on a triangle element in barycentric coordinates.
  QuadPtWt(N_TRI_QUAD,
	   Tri_Quad_Pt,		// output
	   Tri_Quad_Wt);	// output

  // Prathish's code
  Ifpack Factory_VP;
  Ifpack_Preconditioner *Prec_VP;
  string PrecType_VP = "ILU";//Change here for diffrent method (ILU,ILUT,Amesos,LU)
  Prec_VP = Factory_VP.Create(PrecType_VP, &A_VP);
  assert (Prec_VP != 0);
  
  Ifpack Factory_Bio;
  Ifpack_Preconditioner *Prec_Bio;
  string PrecType_Bio = "ILU";
  Prec_Bio = Factory_Bio.Create(PrecType_Bio, &A_Bio);
  assert (Prec_Bio != 0);

  Teuchos::ParameterList List_VP;
  //List_VP.set("fact: level-of-fill", 3);//not for ILU
  //List_VP.set("relaxation: type", "Gauss-Seidel");
  //List_VP.set("fact: drop tolerance", 1.0e-4);
  //List_VP.set("fact: relax value", 0.0);
  //List_VP.set("fact: absolute threshold", 0.0);
  //List_VP.set("fact: relative threshold", 1.0);
  //List_VP.set("amesos: solver type", "Amesos_Lapack");

  Teuchos::ParameterList List_Bio;
  //List_Bio.set("fact: level-of-fill", 3);//not for ILU
  //List_Bio.set("relaxation: type", "Gauss-Seidel");
  //List_Bio.set("fact: drop tolerance", 1.0e-4);
  //List_Bio.set("fact: relax value", 0.0);
  //List_Bio.set("fact: absolute threshold", 0.0);
  //List_Bio.set("fact: relative threshold", 1.0);
  //List_Bio.set("amesos: solver type", "Amesos_Lapack");

  Epetra_LinearProblem LP_VP(&A_VP,&x_VP,&b_VP);

  AztecOO Solver_VP(LP_VP);

  Epetra_LinearProblem LP_Bio(&A_Bio,&x_Bio,&b_Bio);
  
  AztecOO Solver_Bio(LP_Bio);

  Prec_VP->SetParameters(List_VP);
  Prec_VP->Initialize();
  
  Prec_Bio->SetParameters(List_Bio);
  Prec_Bio->Initialize();
  
  Solver_VP.SetPrecOperator(Prec_VP);
  
  Solver_VP.SetAztecOption(AZ_solver,AZ_gmres);
  Solver_VP.SetAztecOption(AZ_kspace, 1000);
  Solver_VP.SetAztecOption(AZ_conv, AZ_noscaled);
  Solver_VP.SetAztecOption(AZ_output,100);
//   Solver_VP.SetAztecOption(AZ_pre_calc, AZ_reuse);
  
  Solver_Bio.SetPrecOperator(Prec_Bio);
  
  Solver_Bio.SetAztecOption(AZ_solver,AZ_gmres);
  Solver_Bio.SetAztecOption(AZ_kspace, 1000);
  Solver_Bio.SetAztecOption(AZ_conv, AZ_noscaled);
  Solver_Bio.SetAztecOption(AZ_output,100);

  // Initialize Diffusion counter
  int Diff_Counter = 3;
  
  for(int Count = 0; Count < Diff_Counter; Count++)
  {
    
    if(myid == 0)
    {
      std::cout << "Diffusion step n = " << Count << endl;
    }

    A_Bio.PutScalar(0.0);
    b_Bio.PutScalar(0.0);
    x_Bio.PutScalar(0.0);
    Bio_Soln_Next_t.PutScalar(0.0);
  
    BioDiffSparseAssembly(BIO_DIFF_COEFF,
			  Time_Step,
			  NX,
			  NY,
			  NGP,
			  Vel_Npe,
			  Nem, 
			  N_TRI_QUAD,
			  Vel_Nnm,
			  Vel_Nod, 
			  Bio_Nod_BC,
			  Vel_Glxy,
			  Tri_Quad_Pt, 
			  Tri_Quad_Wt,
			  GAUSPT,
			  GAUSWT,
			  Bio_Soln_Cur_t.Values(),
			  Proc_Node_Part.All_Proc_Nodes_Bio,
			  Proc_Node_Part.My_Proc_Eles,
			  myid,
			  A_Bio, 
			  b_Bio);

    Prec_Bio->Compute();
    
    Solver_Bio.Iterate(1000, DIFF_TOL);
    
    Bio_Soln_Next_t.Import(x_Bio,CompleteSolution_Importer_Bio,Add);

    Bio_Soln_Cur_t = Bio_Soln_Next_t;
    
    sprintf(BioFilename, "./Data/Biofilm/DiffusionOnly/BioFile_00%d.data", Diff_File_No);
  
    WriteBioData(BioFilename,
		  Bio_Soln_Cur_t.Values(), 
		  2 * NX + 1);
    
    Diff_File_No++;

  }
  
  if(myid == 0)
  {
    std::cout << "Diffusion Complete" << endl;
  }
  
  sprintf(BioFilename, "./Data/Biofilm/AdvecDiff/BioFile_00%d.data", Bio_File_No);
  
  WriteBioData(BioFilename,
		Bio_Soln_Cur_t.Values(), 
		2 * NX + 1);
  
  Bio_File_No++;

  // Initialize time counter and EquibTol
  n = 0;
  EquibTol = 1.0;
  
  while ((n <= MAX_TIME_STEP_NUM) & (EquibTol >= EQUIB_TOL))
  {

    if(myid == 0)
    {
      std::cout << "Time step n = " << n << endl;
    }

    // Initialize counter
    L = 0;
  
    // Initialize tolerances
    Tol = 1.0;

    Soln_Cur_L = Soln_Cur_t;
    Bio_Soln_Next_t = Bio_Soln_Cur_t;
    Str_Old = Str;

    if(myid == 0)
    {
      std::cout << "Begin Departure Foot" << endl;
    }
    
    GaussDepartureFeet(myid,
		      Nem,
		      N_TRI_QUAD,
		      Vel_Npe,
		      VEL_FLAG,
		      Vel_Nnm,
		      Time_Step,
		      Vel_Glxy,
		      Tri_Quad_Pt,
		      Vel_Nod,
		      Ele_Neigh,
		      Proc_Node_Part.My_Proc_Eles,
		      Soln_Cur_t.Values(),
		      GaussDepartFootx,
		      GaussDepartFooty,
		      GaussDepartElement);
    
    if(myid == 0)
    {
      std::cout << "Begin Stress Departure Foot" << endl;
    }
    
    StrNodeDepartureFeet(myid,
			  Nem,
			  Vel_Npe,
			  Str_Npe,
			  VEL_FLAG,
			  Vel_Nnm,
			  NX,
			  STRESS_FLAG,
			  Time_Step,
			  Vel_Glxy,
			  Vel_Nod,
			  Ele_Neigh,
			  Proc_Node_Part.All_Proc_Nodes_VP,
			  Proc_Node_Part.My_Proc_Eles,
			  Soln_Cur_t.Values(),
			  StrNodeDepartFootx,
			  StrNodeDepartFooty,
			  StrNodeDepartElement);
    
    if(myid == 0)
    {
      std::cout << "Finished Stress Departure Foot" << endl;
    }
    
    while (Tol > TOL) // && L <= 10)
    {
      
      A_VP.PutScalar(0.0);
      b_VP.PutScalar(0.0);
      Soln_Next_L.PutScalar(0.0);
      
      if(myid == 0)
      {
	std::cout << "Begin InitialGuess" << endl;
      }
      
      if((n == 0) & (L == 0))
      {
	InitialGuess(Proc_Node_Part.All_Proc_Nodes_VP.Values(), 
		    Vel_Npe,
		    Vel_Nnm,
		    Proc_Node_Part.My_Proc_Eles.Values(),
		    Vel_Nod,
		    Pre_Nod,
		    Nem,
		    myid,
		    Soln_Cur_t.Values(), 
		    x_VP);
      }
      else
      {
	InitialGuess(Proc_Node_Part.All_Proc_Nodes_VP.Values(), 
		    Vel_Npe,
		    Vel_Nnm,
		    Proc_Node_Part.My_Proc_Eles.Values(),
		    Vel_Nod,
		    Pre_Nod,
		    Nem,
		    myid,
		    Soln_Cur_L.Values(), 
		    x_VP);
      }
      
      if(myid == 0)
      {
	std::cout << "Finished InitialGuess" << endl;
      }
      
      if(myid == 0)
      {
	std::cout << "Vel, Pre & Str Iteration step L = " << L << endl;
      }

      // Assembles the complete matrix for velocity and pressure and its solution vector
      SparseAssembly(SOL_NEW_VIS, 
		      POL_NEW_VIS, 
		      DENSITY_WATER,
		      DENSITY_BIO,
		      T_ZERO,
		      L_ZERO,
		      Time_Step,
		      NX,
		      NY,
		      NGP,
		      Vel_Npe, 
		      Pre_Npe, 
		      Str_Npe,
		      Nem, 
		      N_TRI_QUAD,
		      Vel_Nnm,
		      Pre_Nnm,
		      Str_Nnm,
		      Vel_Nod, 
		      Pre_Nod,
		      Vel_Nod_BC_Hor,
		      Vel_Nod_BC_Ver,
		      Pre_Nod_BC,
		      Vel_Glxy, 
		      Pre_Glxy,
		      VEL_FLAG, 
		      PRE_FLAG, 
		      STRESS_FLAG,
		      Tri_Quad_Pt, 
		      Tri_Quad_Wt,
		      GAUSPT,
		      GAUSWT,
		      Str, 
		      Bio_Soln_Cur_t.Values(),
		      Soln_Cur_t.Values(),
		      Proc_Node_Part.All_Proc_Nodes_VP,
		      Proc_Node_Part.My_Proc_Eles,
		      myid,
		      GaussDepartFootx,
		      GaussDepartFooty,
		      GaussDepartElement,
		      A_VP, 	// output
		      b_VP);	// output

      // Prathish's Code
      Prec_VP->Compute();

      Solver_VP.Iterate(10000, TOL);

      Soln_Next_L.Import(x_VP,CompleteSolution_Importer_VP,Add);

      if(myid == 0)
      {
	std::cout << "Vel, Pre Finished" << endl;
      }
      
      Tol = VectorDiffNorm(Soln_Cur_L.Values(), 
			Soln_Next_L.Values(), 
			2 * Vel_Nnm + Pre_Nnm);

      Soln_Cur_L = Soln_Next_L;

      if(myid == 0)
      {
	std::cout << "Str Begun" << endl;
      }
      
      // Determine the stress tensor.
      Conformation(Time_Step,
		   SOL_NEW_VIS,
		   POL_NEW_VIS,
		   NX,
		   RETARD_TIME,
		   DENSITY_WATER,
		   DENSITY_BIO,
		   T_ZERO,
		   L_ZERO,
		   Vel_Npe, 
		   Pre_Npe,
		   Str_Npe, 
		   Nem, 
		   Vel_Nnm,
		   myid,
		   Vel_Nod, 
		   Pre_Nod, 
		   Str_Nod,
		   Vel_Glxy, 
		   Pre_Glxy, 
		   Str_Glxy,
		   Ele_Neigh,
		   VEL_FLAG, 
		   STRESS_FLAG, 
		   Soln_Cur_L.Values(),
		   Soln_Cur_t.Values(),
		   Bio_Soln_Cur_t.Values(),
		   Str,
		   Str_Old,
		   StrNodeDepartFootx,
		   StrNodeDepartFooty,
		   StrNodeDepartElement,
		   Proc_Node_Part.All_Proc_Nodes_VP,
		   Proc_Node_Part.My_Proc_Eles,
		   Str_New);	//output

      // Update the iteration variables
      Str = Str_New;
//       StressReassemble(Str_Nnm,
// 		       myid,
// 		       Str_New, 
// 		       Proc_Node_Part.All_Proc_Nodes_VP, 
// 		       Str);

      if(myid == 0)
      {
	std::cout << "Str Finished" << endl;
      }
      
      // Iterate counter
      L++;

    } // L while loop

    A_Bio.PutScalar(0.0);
    b_Bio.PutScalar(0.0);
    x_Bio.PutScalar(0.0);
    Bio_Soln_Next_t.PutScalar(0.0);

    if(myid == 0)
    {
      std::cout << "Bio step " << endl;
    }

    BioSparseAssembly(DIFF_COEFF,
		      Time_Step,
		      NX,
		      NY,
		      NGP,
		      Vel_Npe,
		      Nem, 
		      N_TRI_QUAD,
		      Vel_Nnm,
		      Vel_Nod, 
		      Bio_Nod_BC,
		      Vel_Glxy,
		      Tri_Quad_Pt, 
		      Tri_Quad_Wt,
		      GAUSPT,
		      GAUSWT,
		      Soln_Cur_L.Values(),
		      Soln_Cur_t.Values(),
		      Bio_Soln_Cur_t.Values(),
		      Proc_Node_Part.All_Proc_Nodes_Bio,
		      Proc_Node_Part.My_Proc_Eles,
		      myid,
		      A_Bio, 
		      b_Bio);

    Prec_Bio->Compute();
    
    Solver_Bio.Iterate(1000, BIO_TOL);
    
    Bio_Soln_Next_t.Import(x_Bio,CompleteSolution_Importer_Bio,Add);
    
    Bio_Soln_Cur_t = Bio_Soln_Next_t;

    if(myid == 0)
    {
      std::cout << "Bio Finished" << endl;
    }
    
    if(myid == 0)
    {
      std::cout << "Begin Writing to file" << endl;
    }
    
    //If the timestep n just computed is to be printed to the data files
    if((n % Write_Time_Steps_Skipped) == 0)
    {
      if(File_No <= 9)
      {
	sprintf(HorVelFilename, "./Data/Velocity/HorVel/HorVelFile_00%d.data", File_No);
	sprintf(VertVelFilename, "./Data/Velocity/VertVel/VertVelFile_00%d.data", File_No);
	sprintf(PreFilename, "./Data/Pressure/PreFile_00%d.data", File_No);
	
	sprintf(StrXXFilename, "./Data/Stress/XX/StrXXFile_00%d.data", File_No);
	sprintf(StrXYFilename, "./Data/Stress/XY/StrXYFile_00%d.data", File_No);
	sprintf(StrYYFilename, "./Data/Stress/YY/StrYYFile_00%d.data", File_No);
      }
      else if ((File_No >= 10) && (File_No < 100))
      {
	sprintf(HorVelFilename, "./Data/Velocity/HorVel/HorVelFile_0%d.data", File_No);
	sprintf(VertVelFilename, "./Data/Velocity/VertVel/VertVelFile_0%d.data", File_No);
	sprintf(PreFilename, "./Data/Pressure/PreFile_0%d.data", File_No);
	
	sprintf(StrXXFilename, "./Data/Stress/XX/StrXXFile_0%d.data", File_No);
	sprintf(StrXYFilename, "./Data/Stress/XY/StrXYFile_0%d.data", File_No);
	sprintf(StrYYFilename, "./Data/Stress/YY/StrYYFile_0%d.data", File_No);
      }
      else if ((File_No >= 100) && (File_No < 1000))
      {
	sprintf(HorVelFilename, "./Data/Velocity/HorVel/HorVelFile_%d.data", File_No);
	sprintf(VertVelFilename, "./Data/Velocity/VertVel/VertVelFile_%d.data", File_No);
	sprintf(PreFilename, "./Data/Pressure/PreFile_%d.data", File_No);
	
	sprintf(StrXXFilename, "./Data/Stress/XX/StrXXFile_%d.data", File_No);
	sprintf(StrXYFilename, "./Data/Stress/XY/StrXYFile_%d.data", File_No);
	sprintf(StrYYFilename, "./Data/Stress/YY/StrYYFile_%d.data", File_No);
      }
      
      if(Bio_File_No <= 9)
      {
	sprintf(BioFilename, "./Data/Biofilm/AdvecDiff/BioFile_00%d.data", Bio_File_No);
      }
      else if ((Bio_File_No >= 10) && (Bio_File_No < 100))
      {
	sprintf(BioFilename, "./Data/Biofilm/AdvecDiff/BioFile_0%d.data", Bio_File_No);
      }
      else if ((Bio_File_No >= 100) && (Bio_File_No < 1000))
      {
	sprintf(BioFilename, "./Data/Biofilm/AdvecDiff/BioFile_%d.data", Bio_File_No);
      }
      
      WriteVelPreData(HorVelFilename, 
		      VertVelFilename,
		      PreFilename , 
		      Soln_Cur_t.Values(), 
		      NX);

      WriteBioData(BioFilename,
		      Bio_Soln_Cur_t.Values(), 
		      2 * NX + 1);

      WriteStrData(StrXXFilename,
		  StrXYFilename,
		  StrYYFilename,
		  Str, 
		  Str_Write_Size);
      
      File_No++;
      Bio_File_No++;
    }
    
    if(myid == 0)
    {
      std::cout << "Finished Writing to file" << endl;
    }
    
//     A_VP.PutScalar(0.0);
//     b_VP.PutScalar(0.0);
//     Soln_Next_L.PutScalar(0.0);
    
//     InitialGuess(Proc_Node_Part.All_Proc_Nodes_VP.Values(), 
// 		  Vel_Npe,
// 		  Vel_Nnm,
// 		  Proc_Node_Part.My_Proc_Eles.Values(),
// 		  Vel_Nod,
// 		  Pre_Nod,
// 		  Nem,
// 		  myid,
// 		  Soln_Cur_t.Values(), 
// 		  x_VP);
    
//     if(myid == 0)
//     {
//       std::cout << "Equib check " << endl;
//     }
//     
//     EquibSparseAssembly(SOL_NEW_VIS, 
// 			POL_NEW_VIS, 
// 			DENSITY_WATER,
// 			DENSITY_BIO,
// 			T_ZERO,
// 			L_ZERO,
// 			Time_Step,
// 			NX,
// 			NY,
// 			NGP,
// 			Vel_Npe, 
// 			Pre_Npe, 
// 			Nem, 
// 			N_TRI_QUAD,
// 			Vel_Nnm,
// 			Pre_Nnm,
// 			Vel_Nod, 
// 			Pre_Nod,
// 			Vel_Nod_BC_Hor,
// 			Vel_Nod_BC_Ver,
// 			Pre_Nod_BC,
// 			Vel_Glxy, 
// 			Pre_Glxy, 
// 			Ele_Neigh,
// 			VEL_FLAG, 
// 			PRE_FLAG, 
// 			Tri_Quad_Pt, 
// 			Tri_Quad_Wt,
// 			GAUSPT,
// 			GAUSWT,
// 			Str, 
// 			Bio_Soln_Cur_t.Values(),
// 			Soln_Cur_t.Values(),
// 			Proc_Node_Part.All_Proc_Nodes_VP,
// 			Proc_Node_Part.My_Proc_Eles,
// 			myid,
// 			A_VP, 	// output
// 			b_VP);
//     
//     Prec_VP->Compute();
// 
//     Solver_VP.Iterate(10000, TOL);
// 
//     Soln_Next_L.Import(x_VP,CompleteSolution_Importer_VP,Add);
// 
//     SolnDiffNorm = VectorDiffNorm(Soln_Cur_t.Values(), 
// 				  Soln_Next_L.Values(), 
// 				  2 * Vel_Nnm + Pre_Nnm);
//     
//     SolnNorm = VectorNorm(Soln_Cur_t.Values(), 
// 			  2 * Vel_Nnm + Pre_Nnm);
//     
    // Dont update Vel and Pre until here so that BioSparseAssembly has the old and new Vels and Pre to use.
    Soln_Cur_t = Soln_Cur_L;
//     
//     EquibTol = SolnDiffNorm / SolnNorm * 1 / (NX * NX);
//     
//     if(myid == 0)
//     {
//       std::cout << "EquibTol = " << EquibTol << endl;
//     }
    
    n++;
    
  } // n while loop

  return 0;
}