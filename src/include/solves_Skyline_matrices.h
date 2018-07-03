
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

  void Skyline(int &NEqM, int &NEl, int &NNode, int &NDOF, int *&NTK, int **&INod, int **&ID, int *&JD);
  void Reduce_Skyline(int &NEqM, double *&K_S, int *&NTK, int *&JD, ofstream &info);
  void Gauss_El_Skyline(int *&NTK, int *&JD, int &NEqM, double *&UN, double *&K_S);
  void Matrix_Multiplication(int *&NTK, int *&JD, double *&M1, double *&M2, double *&M3, int NEqM);

public:
  solve_Skyline_matrices_cls();

}
} // namespace Solver_ns

} // namespace main_ns

#endif
