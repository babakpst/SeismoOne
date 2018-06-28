

#include <iostream>
#include <string>
#include <sstream>

#include "../include/Load.h"
#include "../include/discretize_the_domain_cls.h"
#include "../include/create_global_matrices_cls.h"

#ifndef SOLVER_H
#define SOLVER_H

namespace main_ns
{

namespace Solver_ns
{

class Solver_cls 
{


  main_ns::Solver_ns::apply_seismic_loads_to_the_domain_cls* Loads;
  
protected:
public:
  main_ns::discretization_ns::discretization_cls* DiscretizedModel;
  main_ns::model_ns::model_cls* Model;
  main_ns::Matrices_ns::Matrices_cls* Matrices;

  Solver_cls(main_ns::discretization_ns::discretization_cls*, main_ns::model_ns::model_cls*,
             main_ns::Matrices_ns::Matrices_cls*, 
             main_ns::Solver_ns::apply_seismic_loads_to_the_domain_cls*)

  virtual solve_the_system_using_implicit_newmark_method(void) = 0;

  virtual ~Solver_cls();
};

} // namespace Solver_ns
} // namespace main_ns

#endif