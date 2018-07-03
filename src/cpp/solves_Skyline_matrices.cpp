
#include "../include/solve_Skyline_matrices.h"

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

  for (int i = 0; i < Matrix->JD[DiscretizedModel->NEqM - 1]; i++)
  {
    Matrix->K[i] = Matrix->K[i] + A0 * Matrix->M[i] + A1 * Matrix->C[i];
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

void main_ns::Solver_ns::solve_Skyline_matrices_cls::Matrix_Multiplication(double *&Matrix, double *&Temp, double *&UN)(int *&NTK, int *&JD, double *&M, double *&M2, double *&M3, int NEqM)
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

///////////////////////////////////////////////////////////////////////////////////////////////////
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

void main_ns::Solver_ns::solve_Skyline_matrices_cls::Skyline(int &NEqM, int &NEl, int &NNode, int &NDOF, int *&NTK, int **&INod, int **&ID, int *&JD)
{

  int NEqEl;

  int *ND; // element constraints

  NEqEl = NDOF * NNode;
  ND = new int[NEqEl];

  for (int i = 0; i < NEqEl; i++)
  {
    ND[i] = 0.0;
  }

  for (int i = 0; i < NEqM; i++)
  {
    NTK[i] = i;
  }

  {
    int i, j, k;
    for (int iel = 0; iel < NEl; iel++)
    {
      for (int i = 0; i < NNode; i++)
      {

        k = INod[i][iel];
        for (int j = 0; j < NDOF; j++)
        {
          ND[j * NNode + i] = ID[k][j];
        }
      }

      for (int l = 0; l < NEqEl; l++)
      {
        for (int k = 0; k < NEqEl; k++)
        {
          i = ND[l];
          j = ND[k];
          //if ((i==0) ||  (j==0)) continue;
          if (i > j)
            continue;
          if (i < NTK[j])
            NTK[j] = i;
        }
      }
    }
  }

  JD[0] = 0;

  for (int i = 1; i < NEqM; i++)
  {
    JD[i] = JD[i - 1] + i + 1 - NTK[i];
  }

  std::cout << "End function skyline" << std::endl;
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
void main_ns::Solver_ns::solve_Skyline_matrices_cls::Gauss_El_Skyline(int *&NTK, int *&JD, int &NEqM, double *&UN, double *&K)
{

  N = NEqM;
  N1 = N - 1;

  int k, KJ, N, N1, KK, K1, NN; // temporary variables
  for (int k = 0; k < N1; k++)
  {
    K1 = k + 1;
    for (int j = K1; j < N; j++)
    {
      if (k < NTK[j])
        continue;
      KJ = JD[j] + k - j;
      KK = JD[k];
      UN[j] = UN[j] - UN[k] * K[KJ] / K[KK];
    }
  }

  NN = JD[N];
  UN[N] = UN[N] / K[NN];

  for (int i = 0; i < N1; i++)
  {
    k = N - i;
    K1 = k + 1;
    for (int j = K1; j < N; j++)
    {
      if (k >= NTK[j])
      {
        KJ = JD[j] + k - j;
        UN[k] = UN[k] - K[KJ] * UN[j];
      }
    }
    KK = JD[k];
    UN[k] = UN[k] / K[KK];
  }
}
