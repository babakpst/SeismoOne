
#include "../include/solves_Skyline_matrices.h"

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

void main_ns::Solver_ns::solve_Skyline_matrices_cls::Compute_the_effective_matrix()
{

  std::cout << " Obtaining the effective matrix ..." << std::endl;
  Compute_the_effective_matrix();

  for (int i = 0; i < Matrices->JD[DiscretizedModel->NEqM - 1]; i++)
  {
    Matrices->K[i] = Matrices->K[i] + A0 * Matrices->M[i] + A1 * Matrices->C[i];
  }
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

void main_ns::Solver_ns::solve_Skyline_matrices_cls::Reduce_the_effective_forece()
{

  std::cout << "Reduce effective matrix ..." << std::endl;
  int K1, I1, KJ, KK, KI, IJ; //temp var

  double Fac; //temp var

  int N, N1; // temp var to follow the algorithm

  N = Model->NEqM;
  N1 = N - 1;

  for (int K = 0; K < N1; K++)
  {
    K1 = K + 1;
    for (int J = K1; J <= N; J++)
    {
      if (K < (Matrix->NTK[J]))
        continue;
      KJ = Matrix->JD[J] + K - J;
      KK = Matrix->JD[K];
      I1 = K1;
      if (K1 < (Matrix->NTK[J]))
        I1 = Matrix->NTK[J];
      for (int I = I1; I <= J; I++)
      {
        if (K > (Matrix->NTK[I]))
        {
          KI = Matrix->JD[I] + K - I;
          Fac = Matrix->K[KI] / K[KK];
          IJ = Matrix->JD[J] + I - J;
          Matrix->K[IJ] = Matrix->K[IJ] - Matrix->K[KJ] * Fac;
        }
      }
    }
  }
}

/*
###################################################################################################
Purpose: This function conducts matrix-matrix multiplication.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/28/2018 - Subroutine initiated.
V0.01: 07/02/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/

void main_ns::Solver_ns::solve_Skyline_matrices_cls::
                                 Matrix_Multiplication(double *&Matrix, double *&Temp, double *&UN)
{

  int I, J, K, IJ;

  for (int II = 0; II < DiscretizedModel->NEqM; II++)
  {
    for (int JJ = 0; JJ < DiscretizedModel->NEqM; JJ++)
    {
      I = II;
      J = JJ;
      if (I > J)
      {
        K = I;
        I = J;
        J = K;
      }
      if (I < NTK[J])
        continue;
      IJ = JD[J] + I - J;
      UN[II] += Matrix[IJ] * Temp[JJ];
    }
  }
}

/*
###################################################################################################
Purpose: This function solver the system for each RHS at each time step.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/28/2018 - Subroutine initiated.
V0.01: 07//2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/
void main_ns::Solver_ns::solve_Skyline_matrices_cls::
                              Solve_the_system_for_this_RHS_using_Gaussina_Elimination(double *&UN)
{

  int k, KJ, N, N1, KK, K1, NN; // temporary variables
  
  N = Model->NEqM;
  N1 = N - 1;
 
  for (int k = 0; k < N1; k++)
  {
    K1 = k + 1;
    for (int j = K1; j < N; j++)
    {
      if (k < Matrix->NTK[j])
        continue;
      KJ = Matrix->JD[j] + k - j;
      KK = Matrix->JD[k];
      UN[j] = UN[j] - UN[k] * Matrix->K[KJ] / Matrix->K[KK];
    }
  }

  NN = Matrix->JD[N];
  UN[N] = UN[N] / Matrix->K[NN];

  for (int i = 0; i < N1; i++)
  {
    k = N - i;
    K1 = k + 1;
    for (int j = K1; j < N; j++)
    {
      if (k >= Matrix->NTK[j])
      {
        KJ = Matrix->JD[j] + k - j;
        UN[k] = UN[k] - Matrix->K[KJ] * UN[j];
      }
    }
    KK = Matrix->JD[k];
    UN[k] = UN[k] / Matrix->K[KK];
  }
}
