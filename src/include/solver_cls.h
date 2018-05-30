
/*
Purpose: This class defines all virutal functions related to all solvers.
*/

#include <iostream>
#include <string>
#include <sstream>

#ifndef SOLVER_CLS_H
#define SOLVER_CLS_H


namespace main_ns
{

namespace solver_ns
{

class solver_cls{

  public:
  solver_cls();
  virtual void matrices_fn (void) =0;
  virtual void assemble_fn (void) =0;
  virtual void shapefunctions_fn (void) =0;
  virtual void load_fn (void) =0;
  virtual void solver_fn (void) =0;
  virtual void results_fn (void) =0;



};


}
}
#endif // !SOLVER_CLS_H