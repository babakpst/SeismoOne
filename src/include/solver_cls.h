

#include <iostream>
#include <string>
#include <sstream>

#include "../include/Load.h"

#ifndef SOLVER_H
#define SOLVER_H

namespace main_ns
{

namespace Solver_ns
{

class Solver_cls : public main_ns::Solver_ns::apply_seismic_loads_to_the_domain
{

protected:
public:
  Solver_cls();

  virtual solve_the_system_using_implicit_newmark_method(void) = 0;

  virtual ~Solver_cls();
};

} // namespace Solver_ns
} // namespace main_ns

#endif