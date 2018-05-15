
// include files
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>


using namespace std;


#ifndef VARIABLES_H
#define VARIABLES_H

// = Global Variables ===============================================================================================================================

  // - Integer variables ----------------------------------------------------------------------------------------------------------------------------
extern  int NDim;             // number of Dimension
extern  int NInt;             // number of integration points in Gauss integration scheme
extern  int NNode;            // number of nodes of element
extern  int NDOF;             // number of degrees of freedom
extern  int NMat;             // number of Materials
extern  int NPM;              // number of Properties of Materials
extern  int NStep;            // number of time steps
extern  int NEl;              // number of Elements
extern  int NJ;               // number of Joints(Nodes)
extern  int NEqM;             // number of Equations (Modified i.e. after applying essential boundary conditions)
extern  int OShFunc;          // order of Shape Functions
extern  int LoadType;         // Load Type 0:prssure- 1: DRM
extern  int NNBndry ;         // Number of nodes on the bounday of the DRM
extern  int NNLayer ;         // Number of nodes on the bounday layer of the DRM
extern  int Wave_Type;        // Wave type (0: SV- 1: P) 
extern  int Wave_Func;        // Wave_Func (0: sine-1: Ricker)
extern  int Dis_History;      // Number of nodes for the history of displacement
extern  int Solver;           // Solver type 0: full matrices 1: skyline method
extern  int NNodePWaveL;      // Number of Nodes per wavelength
extern  int NEl_DRM;          // DRM boundary


// - Real variables -------------------------------------------------------------------------------------------------------------------------------
extern  double DT ;           // time step
extern  double Gama, Beta ;   // Newmark parameters
extern  double Alpha, P ;     // loading parameters
extern  double A ;            // area of the cross section of the beam
extern  double L;             // length of the beam
extern  double alpha1, alpha2;// Ricker pulse signal
extern  double amplitude ;    // Amplitude of the incident wave
extern  double omega ;        // central cyclic frequency in the Ricker pulse wave

extern  const double pi;


// - Integer arrays -------------------------------------------------------------------------------------------------------------------------------
extern  int * MTel;           // material type of each Element
extern  int * NoBndry_DRM;    // a vector that holds the node numbers on the DRM boundary
extern  int * NoLayer_DRM;    // a vector that holds the node numbers on the DRM layer
extern  int * JD;             // Skyline matrix
extern  int * NTK;            // Skyline matrix

extern  int * ND_b;           // Nodal ID for DRM
extern  int * ND_e;           // Nodal ID for DRM
 
extern  int * Nodal_History;  // a vector that holds the node numbers to record the history of displacement.
extern  int * Element_Layer;  // Number of elements in each layer
extern  int * Layer_Depth;    // Depth of each layer

extern  int ** INod;          // node connectivity
extern  int ** ID;            // identity

// - Real arrays ----------------------------------------------------------------------------------------------------------------------------------
extern  double * F;           // global force vector
extern  double * Length;      // global force vector
extern  double * K_S;          // global stiffness matrix -skyline
extern  double * C_S;          // global damping matrix -skyline
extern  double * M_S;          // global mass matrix -skyline
extern  double * Loc_History;    // a vector that holds the required locations to store the time history of displacement

extern  double ** XYZ;        // node coordinates
extern  double ** PMat;       // properties of materials
extern  double ** K;          // global stiffness matrix
extern  double ** C;          // global damping matrix
extern  double ** M;          // global mass matrix

extern  double ** K_eb;          // global stiffness matrix
extern  double ** C_eb;          // global damping matrix
extern  double ** M_eb;          // global mass matrix

// - data structure -------------------------------------------------------------------------------------------------------------------------------
// Gauss Points
typedef struct  {
    double XInt[3] ;    // ABSCISSAE
    double WInt[3] ;    // WEIGHTS
  } Gauss ;

// Shape functions
typedef  struct  {     // first order shape functions
    double X1 ;
    double FN[2] ;
  } SF_1D_2N ;

  typedef struct  {     // second order shape functions
    double X1 ;
    double FN[3] ;
  }SF_1D_3N;

// Differentials of shape functions
  typedef struct {     // first order shape functions
    double X1 ;
    double DFXI[2] ;
  }DSF_1D_2N ;

  typedef struct  {     // second order shape functions
    double X1 ;
    double DFXI[3] ;
  } DSF_1D_3N;

#endif
