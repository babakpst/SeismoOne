
#include <cmath>
#include <iostream>
#include <string>

#include "../include/reading_the_model_cls.h"

#ifndef LOAD_H
#define LOAD_H

namespace main_ns
{

namespace Solver_ns
{

class apply_seismic_loads_to_the_domain_cls
{

private:

  struct PointLoad 
  {
  int    Wave_Func;  // determine the shape of the incoming wave
  double amplitude;  // amplitude of the incoming wave
  double Time;       // Time instant that we want to calculate the DRM forces
  double x;          // the location that we want to find the forces
  double c;          // the speed of wave
  double omega;      // the frequency of the incoming wave
  double alpha1;     // the lower bound of the phase of the incoming wave
  double alpha2;     // the upper bound of the phase of the incoming wave
  double u;          // displacement at this node of the DRM
  double v;          // velocity at this node of the DRM
  double a;          // the acceleration at this node of the DRM
  };
 
  PointLoad Load;

  void DRM_PointValues(PointLoad );

public:
  double LoadFactor; // We use this apply pressure loads
  const double pi = {3.14159265359};

  apply_seismic_loads_to_the_domain_cls();
  double LoadFunction(const double, const double, const double);

  void DRM_Loads_Implicit(double &alpha1, double &alpha2, double &Time, int NDim, int NNBndry, int NNLayer, int &Wave_Type, int &Wave_Func, double &amplitude, double &c, double *&UN, double **&XYZ, int *&NoBndry_DRM, int *&NoLayer_DRM, double **&M_eb, double **&C_eb, double **&K_eb, int *&ND_e, int *&ND_b);
  void DRM_PointValues_Freq(double &amplitude, double &x, double &c, double &omega, double &u_R, double &u_I);

  void HistorySolution(int &NJ, double &TIME, double &Alpha, double &P, double &E, double &Rho, double &A, double *&U_EX, double **&XYZ);

} // load class
} // namespace Solver_ns
} // namespace main_ns

#endif
