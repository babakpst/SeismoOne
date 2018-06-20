
/*
Purpose: This class defines all virutal functions related to all solvers.
*/

#include <iostream>
#include <string>
#include <sstream>

#include "../include/Discretization_cls.h"
#include "../include/ShapeFunctions_cls.h"

#include "../include/create_full_matrices_cls.h"
#include "../include/create_skyline_matrices_cls.h"

#include "../include/assemble_local_to_global.h"

#ifndef CREATE_GLOBAL_MATRICES_H
#define CREATE_GLOBAL_MATRICES_H

namespace main_ns
{

namespace Matrices_ns
{

class Matrices_cls
{

  int MType;  // material type
  double E;   // elastic modulus
  double Rho; // density

  void assembling_local_matrices_into_global_matrices_fn();

protected:
  int NEqEl; // Number of equations of each element

  double **K_eb; // global stiffness matrix
  double **C_eb; // global damping matrix
  double **M_eb; // global mass matrix

  double *F; // global force vector

  int LoadFunc; // Load Function  0:DRM

  double **XT; // Array to store coordinates of the element
  double **Ke; // stiffness matrix of each element
  double **Ce; // damping matrix of each element
  double **Me; // mass matrix of each element
  double *Fe;  // element force vector

  void allocating_local_matrices_fn();

  virtual void allocating_global_matrices_fn(void) = 0;

  void assembling_local_matrices_into_global_matrices_fn();

  void compute_elemental_matrices_fn(int, double, double);

public:
  main_ns::discretization_ns::discretization_cls *DiscretizedModel;
  main_ns::model_ns::model_cls *Model;

  main_ns::ShapeFunctions_ns::ShapeFunctions_cls *SF;
  main_ns::Matrices_ns::assemble_local_to_global_cls *assemble;

  Matrices_cls(main_ns::discretization_ns::discretization_cls *, main_ns::model_ns::model_cls *);

  virtual ~Matrices_cls();
};

} // namespace Matrices_ns
} // namespace main_ns
#endif