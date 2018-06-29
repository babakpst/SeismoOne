
#include "../include/solver_cls.h"

#ifndef SOLVER_SKYLINE_SYSTEM_H
#define SOLVER_SKYLINE_SYSTEM_H
namespace main_ns
{

namespace Solver_ns
{

class solve_Skyline_matrices_cls : public main_ns::Solver_ns::Solver_cls
{

  void Skyline(int &NEqM, int &NEl, int &NNode, int &NDOF, int *&NTK, int **&INod, int **&ID, int *&JD);
  void Reduce_Skyline(int &NEqM, double *&K_S, int *&NTK, int *&JD, ofstream &info);
  void Gauss_El_Skyline(int *&NTK, int *&JD, int &NEqM, double *&UN, double *&K_S);
  void Matrix_Multiplication(int *&NTK, int *&JD, double *&M1, double *&M2, double *&M3, int NEqM);

public:
  solve_Skyline_matrices_cls();

  virtual void solve_the_system_using_implicit_newmark_method
  (double &L, int &Wave_Type, int &Wave_Func, int &NStep, int &NEqM, int &LoadType, double &Gama, double &Beta, double &DT, double &Alpha, double *&M_S, double *&C_S, double *&K_S, double *&F, double **&PMat, double **&XYZ, ofstream &FullSol, ofstream &History, int *&ND_e, int *&ND_b, int *&Nodal_History, int *&JD, int *&NTK, ofstream &Check);
}

} // namespace Solver_ns

} // namespace main_ns

#endif
