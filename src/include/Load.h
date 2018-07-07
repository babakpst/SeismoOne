
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
 
// This data type packages all the required data to compute the DRM load at a point
struct InputLoad 
  {
  int NDim;         // The dimension of the model (always 1 in this code)
  int NNBndry;      // Number of nodes on the DRM interface (always 1 in the 1D model)
  int NNLayer;      // Number of nodes on the DRM layer (always 2 in the 1D model)
  int Wave_Type;    // shear or pressure waves, even though in the 1D simulation, this really does not matter
  int Wave_Func;    // determine the shape of the incoming wave

  double x;      // The coordinate
  double u;      // Analytical displacement
  double v;      // Analytical velocity
  double a;      // Analytical acceleration

  double Time;       // Time instant that we want to calculate the DRM forces
  double alpha1;     // the lower bound of the phase of the incoming wave
  double alpha2;     // the upper bound of the phase of the incoming wave
  double amplitude;  // amplitude of the incoming wave
  double c;          // the speed of wave
  double omega;      // the frequency of the incoming wave
  
  int *&NoBndry_DRM;
  int *&NoLayer_DRM;
  int *&ND_e;
  int *&ND_b;

  double *&UN;
  double**& XYZ; 
  double **&M_eb;
  double **&C_eb;
  double **&K_eb;
  };

class apply_seismic_loads_to_the_domain_cls
{

private:

  InputLoad LoadPackage;

  void DRM_PointValues();

public:
  double LoadFactor; // We use this apply pressure loads
  const double pi = {3.14159265359};

  apply_seismic_loads_to_the_domain_cls();
  double LoadFunction(const double, const double, const double);
  void DRM_Loads_Implicit(main_ns::Solver_ns::InputLoad*);
 
  
  

  //void DRM_PointValues_Freq(double &amplitude, double &x, double &c, double &omega, double &u_R, double &u_I);
  //void HistorySolution(int &NJ, double &TIME, double &Alpha, double &P, double &E, double &Rho, double &A, double *&U_EX, double **&XYZ);

}; // load class
} // namespace Solver_ns
} // namespace main_ns

#endif
