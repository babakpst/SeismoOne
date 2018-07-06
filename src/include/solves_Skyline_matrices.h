
#include "../include/solver_cls.h"

#ifndef SOLVER_SKYLINE_SYSTEM_H
#define SOLVER_SKYLINE_SYSTEM_H
namespace main_ns
{

namespace Solver_ns
{

class solve_Skyline_matrices_cls : public main_ns::Solver_ns::Solver_cls
{

  virtual void Compute_the_effective_matrix();
  virtual void Reduce_the_effective_forece();
  virtual void Matrix_Multiplication(double*&Matrix, double*& Temp, double*& UN);
  virtual void Solve_the_system_for_this_RHS_using_Gaussina_Elimination(double*& UN);
  
public:
  solve_Skyline_matrices_cls();
 

};
} // namespace Solver_ns

} // namespace main_ns

#endif
