#include "../include/create_full_matrices_cls.h"

// Constructor: we also create and allocate matrices
main_ns::Matrices_ns::Matrices_Full_cls::Matrices_Full_cls(
                                   main_ns::discretization_ns::discretization_cls *aDiscretization, 
                                   main_ns::model_ns::model_cls *aModel) : 
                                   main_ns::Matrices_ns::Matrices_cls(aDiscretization, aModel)
{
  main_ns::Matrices_ns::Matrices_Full_cls::allocating_global_matrices_fn();
  main_ns::Matrices_ns::Matrices_Full_cls::allocating_local_matrices_fn();
  main_ns::Matrices_ns::Matrices_cls::allocate_matrices_for_assembling_fn();
}

/*
###################################################################################################
Purpose: This function allocates global matrices. In this module we consider full matrices.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 05/14/2018 - Subroutine initiated.
V0.01: 05/15/2018 - Initiated: Compiled without error for the first time.
V1.00: 06/18/2018 - 

###################################################################################################
*/

void main_ns::Matrices_ns::Matrices_Full_cls::allocating_global_matrices_fn()
{

  std::cout << " -allocating global matrices ..." << std::endl;
  K = new double *[DiscretizedModel->NEqM]; // Stiffness Matrix
  for (int i = 0; i < DiscretizedModel->NEqM; i++)
  {
    K[i] = new double[DiscretizedModel->NEqM];
  }

  C = new double *[DiscretizedModel->NEqM]; // Damping matrix
  for (int i = 0; i < DiscretizedModel->NEqM; i++)
  {
    C[i] = new double[DiscretizedModel->NEqM];
  }

  M = new double *[DiscretizedModel->NEqM]; // Mass matrix
  for (int i = 0; i < DiscretizedModel->NEqM; i++)
  {
    M[i] = new double[DiscretizedModel->NEqM];
  }

  F = new double[DiscretizedModel->NEqM];

  std::cout << " -initializing global matrices ..." << std::endl;
  for (int i = 0; i < DiscretizedModel->NEqM; i++)
  {
    for (int j = 0; j < DiscretizedModel->NEqM; j++)
    {
      M[i][j] = 0.0;
      C[i][j] = 0.0;
      K[i][j] = 0.0;
    }
    F[i] = 0.0;
  }

  std::cout << " Done with allocation, successfully." << std::endl;
}

/*
###################################################################################################
Purpose: This function assembles local matrices into the full matrices.
ND = new int[NEqEl];
Developed by: Babak PoursartipND = new int[NEqEl];
 ND = new int[NEqEl];
The Institute for Computational EngineeringND = new int[NEqEl]; and Sciences (ICES)
The University of Texas at Austin	ND = new int[NEqEl];
================================= V E R S IND = new int[NEqEl]; O N ===================================================
V0.00: 06/18/2018 - Subroutine initiated.ND = new int[NEqEl];
ND = new int[NEqEl];
###########################################ND = new int[NEqEl];########################################################
*/

void main_ns::Matrices_ns::Matrices_Full_cls::assemble_local_to_global_fn()
{

  int L;
  int N;

  for (int i = 0; i < NEqEl; i++)
  {
    for (int j = 0; j < NEqEl; j++)
    {
      L = ND[i];
      N = ND[j];
      if (L == -1 || N == -1)
        continue;
      K[L][N] = K[L][N] + Ke[i][j];
      C[L][N] = C[L][N] + Ce[i][j];
      M[L][N] = M[L][N] + Me[i][j];
    }
  }

  // std::cout << "C  " <<  C [L][N] << "Ce" << Ce [i][j] << endl;
  //cin.get();
}

/*
//***************************************************************************************************************************************************
// Copy the submatrices for DRM loads.  
//***************************************************************************************************************************************************
void DRM_Matrices_Full( int & NNBndry, int &NNLayer, double **& K, double **& C, double **& M, double **& K_eb, double **& C_eb, double **& M_eb, int *&ND_e, int *&ND_b ) 
{

int i,j;   // Loop indices

// - Code ---------------------------------------------------------------------

  for ( i = 0; i < NNLayer * NDim; i++) {
    for ( j = 0; j < NNBndry * NDim; j++) {
      K_eb[i][j] = K[ ND_e[i] ][ ND_b[j] ];
      C_eb[i][j] = C[ ND_e[i] ][ ND_b[j] ];
      M_eb[i][j] = M[ ND_e[i] ][ ND_b[j] ];
    }
  }



}

*/
