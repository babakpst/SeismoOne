

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


  void Gaussian(int &NEqM, double *&UN, double **&K);
  void LDLT(int &NEqM, double **&K);
  void Substitute(int &NEqM, double *&UN, double **&K);

public:
  solve_full_matrices_cls();

}

} // namespace Solver_ns

} // namespace main_ns

//void Reduce_Full (int& NEqM, double **& K, ofstream& Check);


#endif
