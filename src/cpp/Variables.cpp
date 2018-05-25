
#include "../include/Variables.h"

// = Global Variables ===============================================================================================================================

  // - Integer variables ----------------------------------------------------------------------------------------------------------------------------
  
  int LoadFunc;         // Load Function  0:DRM


  const double pi = 3.141592653589793 ;

  // - Integer arrays -------------------------------------------------------------------------------------------------------------------------------
  int * JD;             // Skyline matrix
  int * NTK;            // Skyline matrix

  int * ND_b;           // Nodal ID for DRM
  int * ND_e;           // Nodal ID for DRM

  

  // - Real arrays ----------------------------------------------------------------------------------------------------------------------------------

  double * F;           // global force vector

  double * K_S;         // global stiffness matrix -skyline
  double * C_S;         // global damping matrix -skyline
  double * M_S;         // global mass matrix -skyline

  

  double ** K;          // global stiffness matrix
  double ** C;          // global damping matrix
  double ** M;          // global mass matrix

  double ** K_eb;          // global stiffness matrix
  double ** C_eb;          // global damping matrix
  double ** M_eb;          // global mass matrix









