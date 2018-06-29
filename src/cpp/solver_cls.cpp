
#include "../include/solver_cls.h"

main_ns::Solver_ns::Solver_cls::Solver_cls(
    main_ns::discretization_ns::discretization_cls *aDiscretization,
    main_ns::model_ns::model_cls *aModel,
    main_ns::Matrices_ns::Matrices_cls *aMatrices)
    : DiscretizedModel(aDiscretization), Model(aModel), Matrices(aMatrices)
{

  // defining the Newmark constants
  A0 = 1.0 / (Model->Beta * Model->DT * Model->DT);
  A1 = Model->Gama / (Model->Beta * Model->DT);
  A2 = 1.0 / (Model->Beta * Model->DT);
  A3 = 1.0 / (2.0 * Model->Beta) - 1.0;
  A4 = Model->Gama / Model->Beta - 1.0;
  A5 = Model->DT * (Model->Gama / (2.0 * Model->Beta) - 1.0);

  // defining the material properties for the DRM
  E   = Model->PMat[0][0];    // Elastic Modulus of the base material required for the DRM loads
  Rho = Model->PMat[0][1];  // density of the base material required for the DRM loads
  c   = sqrt(E / Rho); // wave velocity of the base material required for the DRM loads

  double *UN;   // temporay arrays for the Newmark algorithm
  double *U;    // temporay arrays for the Newmark algorithm
  double *UD;   // temporay arrays for the Newmark algorithm
  double *UDD;  // temporay arrays for the Newmark algorithm
  double *Temp; // temporay arrays for the Newmark algorithm

  UN   = new double[DiscretizedModel->NEqM];
  U    = new double[DiscretizedModel->NEqM];
  UD   = new double[DiscretizedModel->NEqM];
  UDD  = new double[DiscretizedModel->NEqM];
  Temp = new double[DiscretizedModel->NEqM];

  // Initializing displacement, velocity and acceleration
  for (int i = 0; i < NEqM; i++)
  {
    UN[i]  = 0.0;
    U[i]   = 0.0;
    UD[i]  = 0.0;
    UDD[i] = 0.0;
  }


}