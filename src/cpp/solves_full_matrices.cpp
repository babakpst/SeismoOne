

#include "../include/solves_full_matrices.h"

/*
###################################################################################################
Purpose: This function computes the effective matrix based on the Newmark algorithm.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/29/2018 - Subroutine initiated.
V0.01: 06/29/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/

void main_ns::Solver_ns::solve_full_matrices_cls::Compute_the_effective_matrix()
{

  // Effective stiffness matrix
  std::cout << " Obtaining the effective matrix ..." << std::endl;
  for (int i = 0; i < DiscretizedModel->NEqM; i++)
  {
    for (int j = 0; j < DiscretizedModel->NEqM; j++)
    {
      Matrix->K[i][j] = Matrix->K[i][j] + A0 * Matrix->M[i][j] + A1 * Matrix->C[i][j];
    }
  }
}

/*
###################################################################################################
Purpose: This function reduces the stiffness matrix using the LDLT method.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/28/2018 - Subroutine initiated.
V0.01: 06/29/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/
void main_ns::Solver_ns::solve_full_matrices_cls::Reduce_the_effective_forece()
{

  std::cout << "Reduce effective matrix ..." << std::endl;
  int tempI;
  double *L;

  L = new double[DiscretizedModel->NEqM]; // Identifications

  for (int j = 0; j < DiscretizedModel->NEqM; j++)
  {
    std::cout << j << " reduces out of " << DiscretizedModel->NEqM << std::endl;

    tempI = 5 + j;
    if (tempI > DiscretizedModel->NEqM)
      tempI = DiscretizedModel->NEqM;

    for (int i = j + 1; i < tempI; i++)
    {
      L[i] = Matrix->K[i][j] / Matrix->K[j][j];
    }
    for (int k = j + 1; k < tempI; k++)
    {
      for (int l = j + 1; l < DiscretizedModel->NEqM; l++)
      {
        Matrix->K[k][l] = Matrix->K[k][l] - L[k] * Matrix->K[j][l];
      }
    }
    for (int i = j + 1; i < tempI; i++)
    {
      Matrix->K[i][j] = L[i];
    }
  }
}

/*
###################################################################################################
Purpose: This function solves the AX=B using Gaussian-Elemination method.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/28/2018 - Subroutine initiated.
V0.01: 07/02/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/

void main_ns::Solver_ns::solve_full_matrices_cls::Matrix_Multiplication()
{
  double TempVar;
  for (int i = 0; i < DiscretizedModel->NEqM; i++)
  { // multiplying the vector by the mass matrix
    TempVar = 0.0;
    for (int j = 0; j < DiscretizedModel->NEqM; j++)
    {
      TempVar += Matrix->M[i][j] * Temp[j];
    }
    UN[i] = TempVar;
  }
}




///////////////////////////////////////////
/*
###################################################################################################
Purpose: This function solves the AX=B using Gaussian-Elemination method.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/28/2018 - Subroutine initiated.
V0.01: 06/28/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/
void main_ns::Solver_ns::solve_full_matrices_cls::Gaussian(int &NEqM, double *&UN, double **&K)
{

  int k, l;
  double Fac;
  double sum;

  // Uppertriangular
  std::cout << "Upper" << std::endl;
  for (int j = 0; j < NEqM - 1; j++)
  {
    for (int i = j + 1; i < NEqM - 1; i++)
    {
      Fac = -K[i][j] / K[j][j];
      for (k = j; k < NEqM; k++)
      {
        K[i][k] += Fac * K[j][k];
      }
      UN[i] += Fac * UN[j];
    }
  }

  //Backwardsubstitution
  std::cout << "Backward" << std::endl;
  for (int i = 0; i < NEqM; i++)
  {
    k = NEqM - i - 1;
    sum = 0.0;
    for (int j = 0; j < i; j++)
    {
      l = NEqM - j - 1;
      sum += K[k][l] * UN[l];
    }
    UN[k] = (UN[k] - sum) / K[k][k];
  }
}

/*
###################################################################################################
Purpose: This function solves the AX=B using LDLT method.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/28/2018 - Subroutine initiated.
V0.01: 06/28/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/

void main_ns::Solver_ns::solve_full_matrices_cls::Substitute(int &NEqM, double *&UN, double **&K)
{

  int k, l; // Loop index
  double temp;
  double *L;

  L = new double[NEqM]; // Identifications

  //cout << "Forward" << endl;
  for (int i = 0; i < NEqM; i++)
  {
    temp = 0.0;
    for (int j = 0; j < i; j++)
    {
      temp += K[i][j] * UN[j];
    }
    UN[i] = UN[i] - temp;
  }

  for (int i = 0; i < NEqM; i++)
  {
    UN[i] = UN[i] / K[i][i];
  }

  //cout << "Backward" << endl;
  for (int i = 0; i < NEqM; i++)
  {

    k = NEqM - i - 1;
    temp = 0.0;
    for (int j = 0; j < i; j++)
    {
      l = NEqM - j - 1;
      temp += K[l][k] * L[l];
    }
    L[k] = (UN[k] - temp);
  }

  for (int i = 0; i < NEqM; i++)
  {
    UN[i] = L[i];
  }
}

//// delete
/*
###################################################################################################
Purpose: This function solves the AX=B using LDLT method.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/28/2018 - Subroutine initiated.
V0.01: 06/28/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/
void main_ns::Solver_ns::solve_full_matrices_cls::solve_the_system_using_implicit_newmark_method()
{

  // equivalent to matrix multiplication
  for (int i = 0; i < NEqM; i++)
  {
    TE = 0.0;
    for (int j = 0; j < NEqM; j++)
    {
      TE += C[i][j] * Temp[j];
    }
    UN[i] += TE;
  }

  // Adding load at this time step
  if (LoadType == 0)
  {
    LoadFactor = LoadFunction(Time, Alpha, P); // Pressure load
    for (int ij = 0; ij < NJ; ij++)
    {
      UN[ij] = UN[ij] - F[ij] * LoadFactor;
    }
  }
  else if (LoadType == 1)
  { //        DRM_Load ();
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
  //Gauss_El_Full ( NEqM, UN, K);
  //Gaussian ( NEqM, UN, K);
  Substitute(NEqM, UN, K);

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
