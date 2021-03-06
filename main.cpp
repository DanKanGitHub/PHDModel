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
// #include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_IntSerialDenseMatrix.h"

#include <fstream>
#include <ostream>
#include "VelBioMesh2d.h"
#include "PreMesh2d.h"
#include "StrMesh2d.h"
#include "SparseAssembly.h"
// #include "EquibSparseAssembly.h"
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
#include "StressMatrixDiffNorm.h"
#include "readInBio.h"
#include "readInVelPre.h"
#include "readInStr.h"
#include "WriteParameters.h"

#include "Ifpack_ConfigDefs.h"
// #include "Ifpack_AdditiveSchwarz.h"
// #include "Ifpack_PointRelaxation.h"
// #include "Ifpack_BlockRelaxation.h"
// #include "Ifpack_SparseContainer.h"
// #include "Ifpack_Amesos.h"
// #include "Ifpack_Graph.h"
// #include "Ifpack_Graph_Epetra_CrsGraph.h"
// #include "Ifpack_Graph_Epetra_RowMatrix.h"
// #include "Ifpack_DenseContainer.h"
#include "Ifpack.h"
#include <Ifpack_Preconditioner.h>

typedef Epetra_CrsMatrix E_CM;
typedef Epetra_Map E_Mp;
typedef Epetra_Vector E_V;

typedef Epetra_SerialDenseMatrix E_SDM;
typedef Epetra_IntSerialDenseMatrix E_ISDM;
typedef Epetra_SerialDenseVector E_SDV;
// typedef Epetra_IntSerialDenseVector E_ISDV;

using namespace std;

//Constants
const int NX = 64;			// Number of element intervals in the horizontal direction
const int NY = 64;
const int NGP = 4;			// Numb#include "Epetra_IntSerialDenseVector.h"er of Gauss points in numerical quadrature, used on the boundary
const int N_TRI_QUAD = 7;		// Number of Gauss points in numerical quadrature, used in the element
const int MAX_TIME_STEP_NUM = 1300;	// Maximum number of time interations

const int VEL_FLAG = 2;			// Program flags 1 = linear, 2 = quadratic
const int PRE_FLAG = 1;
const int STRESS_FLAG = 2;

const int Write_Time_Steps_Skipped = 1;

const double XL = 0.0;			// coordinate of left boundary element
const double XR = 1.0;			// coordinate of right boundary element
const double YB = 0.0;			// Coordinate of the bottom boundary of the domain.
const double YT = 1.0;			// Coordinate of the top boundary of the domain.
// Tianyu's Phasse field paper
const double SOL_NEW_VIS = 1.0;	// Solvent Newtonian Viscosity, Water dynamic viscosity at 25 C.
const double POL_NEW_VIS = 1000.0; //10000.0;	// Polymer Newtonian Viscosity, Honey dynamic viscosity at 25 C.

// Issac's Viscoelastic Fluid Description paper
const double RELAX_TIME = 100.0; //100.0 vs 1080;	// Relaxation time, 18min Issac's paper Commonality of Elastic Relaxation Times in Biofilms
						// However, Viscoelastic Fluid Description of Bacterial Bio-film Material Properties has Relax 
						// time on th eorder of 10s and retard time of the order of 100s.
const double RETARD_TIME = 10.0; //10.0 vs 1079;	// Retardation time, !!!must be less than Relaxation time!!!

const double POL_DENSITY = 1.0;	// Density Biofilm, kg/m^3 at 25 C.
const double SOL_DENSITY = 1.0;	// Density Water, kg/m^3 at 25 C.
const double T_SCALE = 60.0;		// seconds (I read Issac's paper for conformation.)
const double L_SCALE = 0.01;		// meters
const double TOL = 0.0000001;		// 0.000001		
const double BIO_TOL = 0.000000000001;
const double DIFF_TOL = 0.000000000001;		// 10^(-12)
// Hits 0.0000001 for 64X64!
const double EQUIB_TOL = 0.0000001;		// The tolerance that determines if the pre-growth simulation has reached equilibrium
const double DIFF_COEFF = 0.001; // 0.01;	// This is a guess.  Diffusion coefffcient for eps diffusion into water
const double BIO_DIFF_COEFF = 0.00001; // 01; // 0.00001
// const double BIO_DIFF_COEFF = 0.001; // 0.00001
const double TIMESTEPSCALAR = 0.05; // 0.01
const double bio_TimeStepScalar = 0.1;
const int timeCounterStartValue = 0; // The last polymer file number
const int useInitalDiffCounter = 7999; // Uses an already diffused polymer

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
  double Tol, Stress_Tol;
  
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
  int NumGlobalElements_VP, NumGlobalElements_Bio;
  double dx, dy, Time_Step, bio_TimeStep, EquibTol, Cur_Time; // SolnNorm, SolnDiffNorm, 
  char HorVelFilename[100];
  char VertVelFilename[100];
  char PreFilename[100];
  char StrXXFilename[100];
  char StrXYFilename[100];
  char StrYYFilename[100];
  char BioFilename[100];
  char ParameterFilename[100];

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
    bio_TimeStep = bio_TimeStepScalar * dx;
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

  // One row for each element and one column for each quadrature point.
  GaussDepartFootx.Shape(Nem, N_TRI_QUAD);
  GaussDepartFooty.Shape(Nem, N_TRI_QUAD);
  GaussDepartElement.Shape(Nem, N_TRI_QUAD);
  
  // One row for each element and one column for each node.
  StrNodeDepartFootx.Shape(Nem, Str_Npe);
  StrNodeDepartFooty.Shape(Nem, Str_Npe);
  StrNodeDepartElement.Shape(Nem, Str_Npe);
  
  // For creating "maps"
  NumGlobalElements_VP = 2 * Vel_Nnm + Pre_Nnm;
  NumGlobalElements_Bio = Vel_Nnm;
  
  // Assigns the correct row length for printing to files depending upon whether
  // stress is approximated using linear or quadratic shape functions, respectively.
//   if(STRESS_FLAG == 1)
//   {
//     Str_Write_Size = NX + 1;
//   }
//   else
//   {
    Str_Write_Size = 2 * NX + 1;
//   }
  
  ProcNodePartitionClass Proc_Node_Part;
  
  // Used to create the "maps" which follow.
  Proc_Node_Part.ProcNodePartition(myid,
				    NumProc, 
				    Vel_Nnm, 
				    Pre_Nnm, 
				    Vel_Npe,
				    Pre_Npe,
				    Vel_Nod,
				    Pre_Nod,
				    Nem);

  // Velocity and Pressure map structures
  E_Mp Local_Proc_Map_VP(-1, Proc_Node_Part.Nodes_Per_Proc_VP, Proc_Node_Part.My_Proc_Nodes_VP.Values(), 0, Comm);
  E_Mp Full_Sol_Map_VP(NumGlobalElements_VP, NumGlobalElements_VP, 0, Comm);
  Epetra_Import CompleteSolution_Importer_VP(Full_Sol_Map_VP, Local_Proc_Map_VP);

  // Polymer map structures
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

  // Initialize the matrix used in the solution A to all zeros for every entry
  // where a value would be incerted into the matrix when it is assembled.
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

  // Output file write counters
  if(readInBio == 0) {
    File_No = 0;
    Diff_File_No = 0;
    Bio_File_No = 0;
  } else {
    File_No = timeCounterStartValue;
    Diff_File_No = timeCounterStartValue;
    Bio_File_No = timeCounterStartValue;
  };
  
  // Parameter Outputs
  if(timeCounterStartValue == 0) {
    sprintf(ParameterFilename, "./Data/ParameterList.txt");
  } else {
    sprintf(ParameterFilename, "./Data/ParameterListForComparison.txt");
  };
  
  WriteParameters(ParameterFilename,
		   NX,
		   NY,
		   NGP,
		   N_TRI_QUAD,
		   MAX_TIME_STEP_NUM,
		   Write_Time_Steps_Skipped,
		   SOL_NEW_VIS,
		   POL_NEW_VIS,
		   RELAX_TIME,
		   RETARD_TIME,
		   POL_DENSITY,
		   SOL_DENSITY,
		   T_SCALE,
		   L_SCALE,
		   DIFF_COEFF,
		   BIO_DIFF_COEFF,
		   TIMESTEPSCALAR);

  // Writes Velocity and Pressure data to files
  if(timeCounterStartValue == 0) {
    
    // Velocity and Pressure file names
    sprintf(HorVelFilename, "./Data/Velocity/HorVel/HorVelFile_00%d.data", File_No);
    sprintf(VertVelFilename, "./Data/Velocity/VertVel/VertVelFile_00%d.data", File_No);
    sprintf(PreFilename, "./Data/Pressure/PreFile_00%d.data", File_No);
    
    WriteVelPreData(HorVelFilename, 
		    VertVelFilename,
		    PreFilename , 
		    Soln_Cur_t.Values(), 
		    NX);
  } else {

    // Velocity and Pressure file names
    sprintf(HorVelFilename, "./Data/Velocity/HorVel/HorVelFile_00%d.data", File_No);
    sprintf(VertVelFilename, "./Data/Velocity/VertVel/VertVelFile_00%d.data", File_No);
    sprintf(PreFilename, "./Data/Pressure/PreFile_00%d.data", File_No);
    
    readInVelPre(NX,
		Soln_Cur_t.Values(), 
		HorVelFilename,
		VertVelFilename,
		PreFilename);
  };

  // Initialize Biofilm vector
  if(timeCounterStartValue == 0) {
    InitialBio(NX,
	      Vel_Glxy,
	      Bio_Soln_Cur_t.Values());
    
    // Polymer File name for the solutions after initial diffusion
//     sprintf(BioFilename, "./Data/Biofilm/AdvecDiff/BioFile_00%d.data", Bio_File_No);
//     
//     // Writes Polymer data to files
//     WriteBioData(BioFilename,
// 		  Bio_Soln_Cur_t.Values(), 
// 		  2 * NX + 1);
    
    // Polymer data filename for initial polymer and the initial diffusion
    sprintf(BioFilename, "./Data/Biofilm/DiffusionOnly/BioFile_00%d.data", Diff_File_No);
    
    // Writes Polymer data to files
    WriteBioData(BioFilename,
		  Bio_Soln_Cur_t.Values(), 
		  2 * NX + 1);
    
  } else {
//     Bio_File_No = timeCounterStartValue; // Temporary assignment for reading
    sprintf(BioFilename, "./Data/Biofilm/AdvecDiff/BioFile_00%d.data", Bio_File_No);

    readInBio(NX,
	      Bio_Soln_Cur_t.Values(), 
	      BioFilename);
  }
  
  if(timeCounterStartValue == 0) {
    // Stress file names
    sprintf(StrXXFilename, "./Data/Stress/XX/StrXXFile_00%d.data", File_No);
    sprintf(StrXYFilename, "./Data/Stress/XY/StrXYFile_00%d.data", File_No);
    sprintf(StrYYFilename, "./Data/Stress/YY/StrYYFile_00%d.data", File_No);
    
    // Writes Stress data to file
    WriteStrData(StrXXFilename,
		StrXYFilename,
		StrYYFilename,
		Str, 
		Str_Write_Size);
  } else {
    
    // Stress file names
    sprintf(StrXXFilename, "./Data/Stress/XX/StrXXFile_00%d.data", File_No);
    sprintf(StrXYFilename, "./Data/Stress/XY/StrXYFile_00%d.data", File_No);
    sprintf(StrYYFilename, "./Data/Stress/YY/StrYYFile_00%d.data", File_No);
    
    readInStr(NX,
	      Str, 
	      StrXXFilename,
	      StrXYFilename,
	      StrYYFilename);
    
  };
  
  if(useInitalDiffCounter != 0) {
    sprintf(BioFilename, "./Data/Biofilm/DiffusionOnly/BioFile_00%d.data", useInitalDiffCounter);

    readInBio(NX,
	      Bio_Soln_Cur_t.Values(), 
	      BioFilename);
  };
  
  // Increment output file counters
  File_No++;
  Diff_File_No++;
  Bio_File_No++;

  // There are at most 3 neightbors in my triangular mesh.
  Ele_Neigh.Shape(Nem,3);

  // Determine which global elements share each side of each element
  // This is not used anymore
//   ElementNeigh(Nem,
// 	       NX,
// 	       Ele_Neigh);	// output

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
  
  Ifpack Factory_Diff;
  Ifpack_Preconditioner *Prec_Diff;
  string PrecType_Diff = "ILU";
  Prec_Diff = Factory_Diff.Create(PrecType_Diff, &A_Bio);
  assert (Prec_Diff != 0);

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
  
  Teuchos::ParameterList List_Diff;

  Epetra_LinearProblem LP_VP(&A_VP,&x_VP,&b_VP);

  AztecOO Solver_VP(LP_VP);

  Epetra_LinearProblem LP_Bio(&A_Bio,&x_Bio,&b_Bio);
  
  AztecOO Solver_Bio(LP_Bio);
  
  Epetra_LinearProblem LP_Diff(&A_Bio,&x_Bio,&b_Bio);
  
  AztecOO Solver_Diff(LP_Diff);

  Prec_VP->SetParameters(List_VP);
  Prec_VP->Initialize();
  
  Prec_Bio->SetParameters(List_Bio);
  Prec_Bio->Initialize();
  
  Prec_Diff->SetParameters(List_Diff);
  Prec_Diff->Initialize();
  
  Solver_VP.SetPrecOperator(Prec_VP);
  
  Solver_VP.SetAztecOption(AZ_solver,AZ_gmres);
  Solver_Bio.SetAztecOption(AZ_precond,AZ_none);
  Solver_VP.SetAztecOption(AZ_kspace, 1000);
  Solver_VP.SetAztecOption(AZ_conv, AZ_noscaled);
  Solver_VP.SetAztecOption(AZ_output,100);
//   Solver_VP.SetAztecOption(AZ_pre_calc, AZ_reuse);
  
  Solver_Bio.SetPrecOperator(Prec_Bio);
  
  Solver_Bio.SetAztecOption(AZ_solver,AZ_gmres);
  Solver_Bio.SetAztecOption(AZ_precond,AZ_none);
  Solver_Bio.SetAztecOption(AZ_kspace, 1000);
  Solver_Bio.SetAztecOption(AZ_conv, AZ_noscaled);
  Solver_Bio.SetAztecOption(AZ_output,100);
  
  Solver_Diff.SetPrecOperator(Prec_Diff);
  
  Solver_Diff.SetAztecOption(AZ_solver,AZ_cg);
//   Solver_Bio.SetAztecOption(AZ_orthog, AZ_classic);
//   Solver_Bio.SetAztecOption(AZ_precond,AZ_sym_GS); // AZ_sym_GS caused a trivial solution, i.e., 0 everywhere
//   Solver_Bio.SetAztecOption(AZ_poly_ord,3);
  Solver_Diff.SetAztecOption(AZ_kspace, 1000);
  Solver_Diff.SetAztecOption(AZ_conv, AZ_noscaled);
  Solver_Diff.SetAztecOption(AZ_output,100);

  // Initialize Diffusion counter
  int Diff_Counter = 8000; // 3;
  
  // "softens" the initial polymer so that the corners to not
  // cause numeric issues
  for(int Count = 0; Count < Diff_Counter; Count++)
  {
    
    if(myid == 0)
    {
      std::cout << "Diffusion step n = " << Count << endl;
    }

    A_Bio.PutScalar(0.0);
    b_Bio.PutScalar(0.0);
    x_Bio.PutScalar(0.0);
    Soln_Cur_t.PutScalar(0.0);
    Bio_Soln_Next_t.PutScalar(0.0);
    
//     if(Count > 100 ) {
//       BIO_DIFF_COEFF = 0.00000001;
//     }
  
    BioDiffSparseAssembly(BIO_DIFF_COEFF,
			  bio_TimeStep,
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

    Prec_Diff->Compute();
    
    Solver_Diff.Iterate(1000, DIFF_TOL);
    
    Bio_Soln_Next_t.Import(x_Bio,CompleteSolution_Importer_Bio,Add);

    Bio_Soln_Cur_t = Bio_Soln_Next_t;
    
//     sprintf(BioFilename, "./Data/Biofilm/DiffusionOnly/BioFile_00%d.data", Diff_File_No);
    
    // Only prints to file every 50 time steps
    if((Count % 500) == 0)
    {
      if(Diff_File_No <= 9)
      {
	sprintf(BioFilename, "./Data/Biofilm/DiffusionOnly/BioFile_00%d.data", Diff_File_No);
      }
      else if ((Diff_File_No >= 10) && (Diff_File_No < 100))
      {
	sprintf(BioFilename, "./Data/Biofilm/DiffusionOnly/BioFile_0%d.data", Diff_File_No);
      }
      else if ((Diff_File_No >= 100) && (Diff_File_No < 1000))
      {
	sprintf(BioFilename, "./Data/Biofilm/DiffusionOnly/BioFile_%d.data", Diff_File_No);
      }
    
      WriteBioData(BioFilename,
		    Bio_Soln_Cur_t.Values(), 
		    2 * NX + 1);
      
      Diff_File_No++;
    }

  }

  if(myid == 0)
  {
    std::cout << "Diffusion Complete" << endl;
  }
  
  if(Diff_File_No <= 9)
  {
    sprintf(BioFilename, "./Data/Biofilm/AdvecDiff/BioFile_00%d.data", Bio_File_No);
  }
  else if ((Diff_File_No >= 10) && (Diff_File_No < 100))
  {
    sprintf(BioFilename, "./Data/Biofilm/AdvecDiff/BioFile_0%d.data", Bio_File_No);
  }
  else if ((Diff_File_No >= 100) && (Diff_File_No < 1000))
  {
    sprintf(BioFilename, "./Data/Biofilm/AdvecDiff/BioFile_%d.data", Bio_File_No);
  }

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
    Stress_Tol = 1.0;

    Soln_Cur_L = Soln_Cur_t;
    Bio_Soln_Next_t = Bio_Soln_Cur_t;
    Str_Old = Str;
    
    // Update to the current time
    Cur_Time = Time_Step * (n + 1);

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
		      NX,
		      Time_Step,
		      dx,
		      dy,
		      Vel_Glxy,
		      Tri_Quad_Pt,
		      Vel_Nod,
		      Ele_Neigh,
		      Proc_Node_Part.My_Proc_Eles,
		      Soln_Cur_t.Values(),
		      GaussDepartFootx,
		      GaussDepartFooty,
		      GaussDepartElement);
    
//     std::cout << "GaussDepartFootx = " << GaussDepartFootx << endl;
//     std::cout << "GaussDepartFooty = " << GaussDepartFooty << endl;
//     std::cout << "GaussDepartElement = " << GaussDepartElement << endl;
    
//     int QWERY;
//     std::cin >> QWERY;
    
    if(myid == 0)
    {
      std::cout << "Begin Stress Departure Foot" << endl;
//       std::cout << "Ele_Neigh = " << Ele_Neigh  << endl;
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
			  dx,
			  dy,
			  Vel_Glxy,
			  Vel_Nod,
			  Ele_Neigh,
			  Proc_Node_Part.All_Proc_Nodes_VP,
			  Proc_Node_Part.My_Proc_Eles,
			  Soln_Cur_t.Values(),
			  StrNodeDepartFootx,
			  StrNodeDepartFooty,
			  StrNodeDepartElement,
			  Str_Nod,
			  Str_Glxy
// 			  Str, 
// 			  Str_Old
			  );
    
//     std::cout << "StrNodeDepartFootx = " << StrNodeDepartFootx << endl;
//     std::cout << "StrNodeDepartFooty = " << StrNodeDepartFooty << endl;
//     std::cout << "StrNodeDepartElement = " << StrNodeDepartElement << endl;
    
//     int QWERY;
//     std::cin >> QWERY;
    
    if(myid == 0)
    {
      std::cout << "Finished Stress Departure Foot" << endl;
    }
    
//     if(n == 10) {
//       Time_Step = 0.1 * dx;
//     }
    
    while ((Tol > TOL) & (Stress_Tol > TOL)) // && L <= 10)
    {
      
      A_VP.PutScalar(0.0);
      b_VP.PutScalar(0.0);
      Soln_Next_L.PutScalar(0.0);
      
      if(myid == 0)
      {
	std::cout << "Begin InitialGuess" << endl;
      }
      
      if(L == 0)
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
		      SOL_DENSITY,
		      POL_DENSITY,
		      T_SCALE,
		      L_SCALE,
		      Time_Step, // Time_Step,
		      Cur_Time, // For initial conditions
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

      if(myid == 0)
      {
	std::cout << "Out of Sparse" << endl;
      }
      
      // Prathish's Code
      Prec_VP->Compute();

      Solver_VP.Iterate(10000, TOL);

      Soln_Next_L.Import(x_VP,CompleteSolution_Importer_VP,Add);

//       Soln_Next_L = Soln_Cur_L;

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
		   SOL_DENSITY,
		   POL_DENSITY,
		   T_SCALE,
		   L_SCALE,
		   Vel_Npe, 
		   Str_Npe, 
		   Nem, 
		   Vel_Nnm,
		   myid,
		   Vel_Nod, 
		   Str_Nod,
		   Vel_Glxy, 
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

//       Stress_Tol = StressMatrixDiffNorm(Str_Nnm,
// 		      myid,
// 		      Str_New,
// 		      Str);
      
      Stress_Tol = 0.0;

      for(int i = 0; i < Str_Nnm; i++)
      {
	Stress_Tol += (Str(i,0) - Str_New(i,0)) * (Str(i,0) - Str_New(i,0));
	Stress_Tol += (Str(i,1) - Str_New(i,1)) * (Str(i,1) - Str_New(i,1));
	Stress_Tol += (Str(i,2) - Str_New(i,2)) * (Str(i,2) - Str_New(i,2));
      }
      
      Stress_Tol = sqrt(Stress_Tol);
      
//       std::cout << "Stress_Tol = " << Stress_Tol << endl;
      
//       Update the iteration variables
      Str = Str_New;
      StressReassemble(Str_Nnm,
		       myid,
		       Str_New, 
		       Proc_Node_Part.All_Proc_Nodes_VP, 
		       Str);

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
// 			SOL_DENSITY,
// 			POL_DENSITY,
// 			T_SCALE,
// 			L_SCALE,
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