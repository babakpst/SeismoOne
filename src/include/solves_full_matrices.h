

#include "../include/solver_cls.h"

#ifndef SOLVER_FULL_SYSTEM_H
#define SOLVER_FULL_SYSTEM_H

namespace main_ns
{

namespace Solver_ns
{

class solve_full_matrices : public main_ns::Solver_ns::Solver_cls
{

  public:
  solve_full_matrices();
  
  virtual void solve_the_system_using_implicit_newmark_method();
}

} // namespace Solver_ns

} // namespace main_ns

//void Reduce_Full (int& NEqM, double **& K, ofstream& Check);
void LDLT(int &NEqM, double **&K);
void Substitute(int &NEqM, double *&UN, double **&K);
//void Gaussian ( int& NEqM, double *& UN, double **& K);
void Newmark_Full(double &L, int &Wave_Type, int &Wave_Func, int &NStep, int &NEqM, int &LoadType, double &Gama, double &Beta, double &DT, double &Alpha, double **&M, double **&C, double **&K, double *&F, double **&PMat, double **&XYZ, ofstream &FullSol, ofstream &History, int *&ND_e, int *&ND_b, int *&Nodal_History);

#endif
