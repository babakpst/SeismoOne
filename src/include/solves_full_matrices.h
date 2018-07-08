

#include "../include/solver_cls.h"

#ifndef SOLVER_FULL_SYSTEM_H
#define SOLVER_FULL_SYSTEM_H

namespace main_ns
{
namespace Solver_ns
{

class solve_full_matrices_cls : public main_ns::Solver_ns::Solver_cls
{

  virtual void Compute_the_effective_matrix();
  virtual void Reduce_the_effective_forece();
  virtual void Matrix_Multiplication(double **&Matrix, double *&Temp, double *&UN);
  virtual void Solve_the_system_for_this_RHS_using_Gaussina_Elimination(double *&UN);

public:
  main_ns::Matrices_ns::Matrices_Full_cls *Matrices;

  solve_full_matrices_cls(main_ns::address_ns::address_cls *, main_ns::model_ns::model_cls *,
                          main_ns::discretization_ns::discretization_cls *,
                          main_ns::Matrices_ns::Matrices_Full_cls *);
};

} // namespace Solver_ns

} // namespace main_ns

//void Reduce_Full (int& NEqM, double **& K, ofstream& Check);

#endif
