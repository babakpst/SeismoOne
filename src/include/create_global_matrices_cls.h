
/*
Purpose: This class defines all virutal functions related to all solvers.
*/

#include <iostream>
#include <string>
#include <sstream>

#include "../include/Discretization_cls.h"
#include "../include/ShapeFunctions_cls.h"

#include "../include/assemble_local_to_global.h"

#ifndef SOLVER_CLS_H
#define SOLVER_CLS_H

namespace main_ns
{

namespace Matrices_ns
{

class Matrices_cls{
  protected:


    int NEqEl;              // Number of equations of each element

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



    int LoadFunc;         // Load Function  0:DRM

    double ** XT;           // Array to store coordinates of the element
    double ** Ke;           // stiffness matrix of each element
    double ** Ce;           // damping matrix of each element
    double ** Me;           // mass matrix of each element
    double * Fe;            // element force vector



    virtual void allocating_global_matrices_fn (void) =0;   
    virtual void allocating_local_matrices_fn(void) = 0;
    virtual void assembling_local_matrices_into_global_matrices_fn(void) = 0;  

    void compute_elemental_matrices_fn(int, double, double);

  public:
    main_ns::discretization_ns::discretization_cls* DiscretizedModel;
    main_ns::model_ns::model_cls* Model;

    main_ns::ShapeFunctions_ns::ShapeFunctions_cls* SF;

    Matrices_cls( main_ns::discretization_ns::discretization_cls*,main_ns::model_ns::model_cls*);
    
    virtual ~Matrices_cls();
  


  
  //virtual void matrices_fn (void) =0;
  //virtual void compute_elemental_matrices_fn (void) =0;
  //virtual void shapefunctions_fn (void) =0;
  //virtual void load_fn (void) =0;
  //virtual void solver_fn (void) =0;
  //virtual void results_fn (void) =0;



};


}
}
#endif