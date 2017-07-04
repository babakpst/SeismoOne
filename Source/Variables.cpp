
#include "Variables.h"

// = Global Variables ===============================================================================================================================

  // - Integer variables ----------------------------------------------------------------------------------------------------------------------------
  int NDim;             // number of Dimension
  int NInt;             // number of integration points in Gauss integration scheme
  int NNode;            // number of nodes of element
  int NDOF;             // number of degrees of freedom
  int NMat;             // number of Materials
  int NPM;              // number of Properties of Materials
  int NStep;            // number of time steps
  int NEl;              // number of Elements
  int NJ;               // number of Joints(Nodes)
  int NEqM;             // number of Equations (Modified i.e. after applying essential boundary conditions)
  int OShFunc;          // order of Shape Functions
  int LoadType;         // Load Type 0:prssure- 1: DRM
  int LoadFunc;         // Load Function  0:DRM
  int NNBndry = 1;      // Number of nodes on the bounday of the DRM
  int NNLayer = 2;      // Number of nodes on the bounday layer of the DRM
  int Wave_Type;        // Wave type (0: SV- 1: P) 
  int Wave_Func;        // Wave_Func (0: sine-1: Ricker)
  int Dis_History;      // Number of nodes for the history of displacement
  int Solver;           // Solver type 0: full matrices 1: skyline method
  int NNodePWaveL;      // Number of Nodes per wavelength
  int NEl_DRM;          // DRM boundary

  // - Real variables -------------------------------------------------------------------------------------------------------------------------------
  double DT ;           // time step
  double Gama, Beta ;   // Newmark parameters
  double Alpha, P ;     // loading parameters
  double A ;            // area of the cross section of the beam
  double L;             // length of the beam
  double alpha1, alpha2;// 
  double amplitude ;    // Amplitude of the incident wave
  double omega ;        // central cyclic frequency in the Ricker pulse wave

  const double pi = 3.141592653589793 ;

  // - Strings --------------------------------------------------------------------------------------------------------------------------------------

  string TempS;         // temporary variable for reading strings from input files
  string Name;          // name of the input file
  string Directory;     // Input/output directory

  string Input_Dir;            // Input directory
  string OutputMatlab_Dir;     // The directory to write down the input file for Matlab visualizer interface
  string Info_Dir;             // The directory to write down general information
  string FullFile_Dir;         // The directory to write down the full results in the time domain analysis
  string HistoryFile_Dir;      // The directory to write down the time history of displacement in the time domain analysis
  string TransferFunction_Dir; // The directory to write down the frequency domain results

  // - bool -----------------------------------------------------------------------------------------------------------------------------------------

  // - Integer arrays -------------------------------------------------------------------------------------------------------------------------------
  int * MTel;           // material type of each Element
  int * NoBndry_DRM;    // a vector that holds the node numbers on the DRM boundary
  int * NoLayer_DRM;    // a vector that holds the node numbers on the DRM layer
  int * JD;             // Skyline matrix
  int * NTK;            // Skyline matrix

  int * ND_b;           // Nodal ID for DRM
  int * ND_e;           // Nodal ID for DRM

  int * Nodal_History;  // a vector that holds the node numbers to record the history of displacement. 
  int * Element_Layer;  // Number of elements in each layer
  int * Layer_Depth;    // Depth of each layer

  int ** INod;          // node connectivity
  int ** ID;            // identity

  // - Real arrays ----------------------------------------------------------------------------------------------------------------------------------
  double * F;           // global force vector
  double * Length;      // global force vector
  double * K_S;         // global stiffness matrix -skyline
  double * C_S;         // global damping matrix -skyline
  double * M_S;         // global mass matrix -skyline
  double * Loc_History; // a vector that holds the required locations to store the time history of displacement
  
  double ** XYZ;        // node coordinates
  double ** PMat;       // properties of materials
  double ** K;          // global stiffness matrix
  double ** C;          // global damping matrix
  double ** M;          // global mass matrix

  double ** K_eb;          // global stiffness matrix
  double ** C_eb;          // global damping matrix
  double ** M_eb;          // global mass matrix




