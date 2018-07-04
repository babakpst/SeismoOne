
#include "../include/solver_cls.h"

main_ns::Solver_ns::Solver_cls::Solver_cls(
    main_ns::discretization_ns::discretization_cls *aDiscretization,
    main_ns::model_ns::model_cls *aModel,
    main_ns::Matrices_ns::Matrices_cls *aMatrices)
    : DiscretizedModel(aDiscretization), Model(aModel), Matrices(aMatrices)
{
}

/*
###################################################################################################
Purpose: This function reduces of the effective stiffness matrix to a triangular one for 
         Gauss-Elimination method using the skyline method

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/28/2018 - Subroutine initiated.
V0.01: 06/28/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/
void main_ns::Solver_ns::Solver_cls::solve_the_system_using_implicit_newmark_method()
{

  std::cout << " -Time integration using the Newmark method for the skyline system-" << std::endl;

  // defining the Newmark constants
  A0 = 1.0 / (Model->Beta * Model->DT * Model->DT);
  A1 = Model->Gama / (Model->Beta * Model->DT);
  A2 = 1.0 / (Model->Beta * Model->DT);
  A3 = 1.0 / (2.0 * Model->Beta) - 1.0;
  A4 = Model->Gama / Model->Beta - 1.0;
  A5 = Model->DT * (Model->Gama / (2.0 * Model->Beta) - 1.0);

  // defining the material properties for the DRM
  E = Model->PMat[0][0];   // Elastic Modulus of the base material required for the DRM loads
  Rho = Model->PMat[0][1]; // density of the base material required for the DRM loads
  c = sqrt(E / Rho);       // wave velocity of the base material required for the DRM loads

  UN   = new double[DiscretizedModel->NEqM];
  U    = new double[DiscretizedModel->NEqM];
  UD   = new double[DiscretizedModel->NEqM];
  UDD  = new double[DiscretizedModel->NEqM];
  Temp = new double[DiscretizedModel->NEqM];

  F    = new double [DiscretizedModel->NEqM] ;

  // Initializing displacement, velocity and acceleration
  for (int i = 0; i < DiscretizedModel->NEqM; i++)
  {
    UN[i]  = 0.0;
    U[i]   = 0.0;
    UD[i]  = 0.0;
    UDD[i] = 0.0;
    F[i]   = 0.0;
  }

  // effective force
  Compute_the_effective_matrix();

  // reduce the effective force
  Reduce_the_effective_forece();

  bool InitialTime = false; // var to check the if the initial time is correct.

  Initial_Time = -Model->L / c;
  std::cout << "Initial time = " << Initial_Time << std::endl;
  
  // Solve the PDE for each step 
  for (int IStep = 0; IStep < Model->NStep + 1; IStep++)
  {

    Elapsed_Time = IStep * Model->DT;
    Time = Initial_Time + IStep * Model->DT; // Time STEP
    std::cout << "Time Step:  " << IStep << "  Time: " << Time << "  Total time:" << 
                                                                        Elapsed_Time << std::endl;

    // Update displacements, velocity and acceleration
    for (int i = 0; i < DiscretizedModel->NEqM; i++)
    {
      Temp[i] = UD[i];
      U[i]    = UN[i] - U[i];
      UD[i]   = A1 * U[i] - A4 * Temp[i] - A5 * UDD[i];
      UDD[i]  = A0 * U[i] - A2 * Temp[i] - A3 * UDD[i];
      U[i]    = UN[i];
      UN[i]   = 0.0;
    }

    // Effective force - stored in UN
    for (int i = 0; i < DiscretizedModel->NEqM; i++)
    { // find the coefficient of the M matrix
      Temp[i] = A0 * U[i] + A2 * UD[i] + A3 * UDD[i];
    }


    // up to here
    // Multiply the mass matrix by the load vector
    Matrix_Multiplication(Matrices->M, Temp, UN);

    for (int i = 0; i < DiscretizedModel->NEqM; i++)
    {
      Temp[i] = A1 * U[i] + A4 * UD[i] + A5 * UDD[i];
    }

    Matrix_Multiplication(Matrices->C, Temp, UN);

    // Adding loads at this time step
    if (Model->LoadType == 0) // Pressure load
    {
      LoadFactor = LoadFunction(Time, Model->Alpha, Model->P); 
      for (int ij = 0; ij < Model->NJ; ij++)
      {
        UN[ij] = UN[ij] - F[ij] * LoadFactor;
      }
    }
    else if (Model->LoadType == 1) //        DRM_Load ();
    { 
      F[0] = 0;

      DRM_Loads_Implicit(alpha1, alpha2, Time, NDim, NNBndry, NNLayer, Wave_Type, Wave_Func, amplitude, c, UN, XYZ, NoBndry_DRM, NoLayer_DRM, M_eb, C_eb, K_eb, ND_e, ND_b);
    }











    // Check whether the initial time is small enough
    if (IStep == 0)
    {
      for (int i = 0; i < NEqM; i++)
      {
        if (UN[i] != 0)
        {
          InitialTime = true;
          break;
        }
      }
      if (InitialTime == true)
      {
        std::cout << "WARNING: REDUCE THE INITIAL TIME" << std::endl;
        return;
      }
    }

    // SOLVE
    Gauss_El_Skyline(NTK, JD, NEqM, UN, K);

    // time history of the solution at some particular nodes
    History << setw(6) << Total_Time;
    for (int i = 0; i < Dis_History; i++)
    {
      History << setw(20) << UN[Nodal_History[i]];
    }
    History << std::endl;

    //FullSol << "Time= " << Time << endl;
    for (int i = 0; i < NEqM; i++)
    {
      FullSol << UN[i] << setw(20);
    }
    FullSol << std::endl;
  }
}