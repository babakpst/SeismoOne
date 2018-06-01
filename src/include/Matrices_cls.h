
/*
Purpose: This class defines all virutal functions related to all solvers.
*/

#include <iostream>
#include <string>
#include <sstream>

#include "../include/Discretization_cls.h"

#ifndef SOLVER_CLS_H
#define SOLVER_CLS_H


namespace main_ns
{

namespace Matrices_ns
{

class Matrices_cls{
  protected:
  int * JD;             // Skyline matrix
  int * NTK;            // Skyline matrix

  double * K_S;         // global stiffness matrix -skyline
  double * C_S;         // global damping matrix -skyline
  double * M_S;         // global mass matrix -skyline

  double ** K;          // global stiffness matrix
  double ** C;          // global damping matrix
  double ** M;          // global mass matrix

  double ** K_eb;          // global stiffness matrix
  double ** C_eb;          // global damping matrix
  double ** M_eb;          // global mass matrix

  double * F;           // global force vector

  int * ND_b;           // Nodal ID for DRM
  int * ND_e;           // Nodal ID for DRM

  int LoadFunc;         // Load Function  0:DRM

  public:
  main_ns::discretization_ns::discretization_cls* DiscretizedModel;
  main_ns::model_ns::model_cls* Model;


	explicit Matrices_cls(main_ns::discretization_ns::discretization_cls*,
             main_ns::model_ns::model_cls*);
  
  virtual void allocating_global_matrices_fn (void) =0;
  //virtual void matrices_fn (void) =0;
  //virtual void assemble_fn (void) =0;
  //virtual void shapefunctions_fn (void) =0;
  //virtual void load_fn (void) =0;
  //virtual void solver_fn (void) =0;
  //virtual void results_fn (void) =0;



};


}
}
#endif // !SOLVER_CLS_H