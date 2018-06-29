

#include <iostream>
#include <string>
#include <sstream>

#include "../include/Load.h"
#include "../include/discretize_the_domain_cls.h"
#include "../include/create_global_matrices_cls.h"

#include "../include/create_full_matrices_cls.h"
#include "../include/create_skyline_matrices_cls.h"

#ifndef SOLVER_H
#define SOLVER_H

namespace main_ns
{

namespace Solver_ns
{

class Solver_cls
{

  // Newmark constants
  double A0;
  double A1;
  double A2;
  double A3;
  double A4;
  double A5;

  double E;   // Elastic Modulus of the base material required for the DRM loads
  double Rho; // density of the base material required for the DRM loads
  double c;   // wave velocity of the base material required for the DRM loads

  double Elapsed_Time; // The elapsed time in the simulation
  double Time;         // The actual simulation time, considering the effects of DRM
  double Initial_Time; // Starting time of the simulation, usually negative, because of the DRM 


  double *UN;   // temporay arrays for the Newmark algorithm
  double *U;    // temporay arrays for the Newmark algorithm
  double *UD;   // temporay arrays for the Newmark algorithm
  double *UDD;  // temporay arrays for the Newmark algorithm
  double *Temp; // temporay arrays for the Newmark algorithm

  main_ns::Solver_ns::apply_seismic_loads_to_the_domain_cls *Loads;


  virtual void Compute_the_effective_matrix(void)=0;
  virtual void Reduce_the_effective_forece(void)=0;

protected:
public:
  main_ns::discretization_ns::discretization_cls *DiscretizedModel;
  main_ns::model_ns::model_cls *Model;
  main_ns::Matrices_ns::Matrices_cls *Matrices;

  Solver_cls(main_ns::discretization_ns::discretization_cls *, main_ns::model_ns::model_cls *,
             main_ns::Matrices_ns::Matrices_cls *,
             main_ns::Solver_ns::apply_seismic_loads_to_the_domain_cls *);

  void solve_the_system_using_implicit_newmark_method(); //either using the skyline or full system
  
  //(double &L, int &Wave_Type, int &Wave_Func, int &NStep, int &NEqM, int &LoadType, double &Gama, double &Beta, double &DT, double &Alpha, double **&M, double **&C, double **&K, double *&F, double **&PMat, double **&XYZ, ofstream &FullSol, ofstream &History, int *&ND_e, int *&ND_b, int *&Nodal_History);

  virtual ~Solver_cls();
};

} // namespace Solver_ns
} // namespace main_ns

#endif