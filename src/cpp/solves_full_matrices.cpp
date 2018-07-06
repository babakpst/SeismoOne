

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

void main_ns::Solver_ns::solve_full_matrices_cls::Matrix_Multiplication(double **&Matrix, double *&Temp, double *&UN)
{
  double TempVar;
  for (int i = 0; i < DiscretizedModel->NEqM; i++)
  { // multiplying the vector by the mass matrix
    TempVar = 0.0;
    for (int j = 0; j < DiscretizedModel->NEqM; j++)
    {
      TempVar += Matrix[i][j] * Temp[j];
    }
    UN[i] = TempVar;
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
V0.01: 07/06/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/

void main_ns::Solver_ns::solve_full_matrices_cls::
                              Solve_the_system_for_this_RHS_using_Gaussina_Elimination(double *&UN)
{

  int k, l; // temporary variables
  
  double temp;
  
  double *L;

  L = new double[Model->NEqM]; // Identifications

  //cout << "Forward" << endl;
  for (int i = 0; i < Model->NEqM; i++)
  {
    temp = 0.0;
    for (int j = 0; j < i; j++)
    {
      temp += Matrix->K[i][j] * UN[j];
    }
    UN[i] = UN[i] - temp;
  }

  for (int i = 0; i < Model->NEqM; i++)
  {
    UN[i] = UN[i] / Matrix->K[i][i];
  }

  //cout << "Backward" << endl;
  for (int i = 0; i < Model->NEqM; i++)
  {

    k = Model->NEqM - i - 1;
    temp = 0.0;
    for (int j = 0; j < i; j++)
    {
      l = Model->NEqM - j - 1;
      temp += K[l][k] * L[l];
    }
    L[k] = (UN[k] - temp);
  }

  for (int i = 0; i < Model->NEqM; i++)
  {
    UN[i] = L[i];
  }
}

