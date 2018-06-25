

#include <iostream>
#include <string>
#include <sstream>


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

  Solver_cls();

  virtual solve_the_system_using_implicit_newmark_method(void)=0;

  virtual ~Solver_cls();
};

} // namespace Solver_ns
} // namespace main_ns

#endif