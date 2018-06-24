

#include <iostream>
#include <string>
#include <sstream>


#include "../include/create_global_matrices_cls.h"

#ifndef SOLVER_H
#define SOLVER_H

namespace main_ns
{

namespace Solver_ns
{

class Solver_cls
{


protected:

public:

  Solver_cls(main_ns::discretization_ns::discretization_cls *, main_ns::model_ns::model_cls *);

  virtual ~Solver_cls();
};

} // namespace Solver_ns
} // namespace main_ns

#endif