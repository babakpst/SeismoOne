#include "../include/create_skyline_matrices_cls.h"

// Constructor: we also create and allocate matrices
main_ns::Matrices_ns::Matrices_Skyline_cls::Matrices_Skyline_cls(main_ns::discretization_ns::discretization_cls *aDiscretization,
                                                                 main_ns::model_ns::model_cls *aModel) : main_ns::Matrices_ns::Matrices_cls(aDiscretization, aModel)
{
  main_ns::Matrices_ns::Matrices_Skyline_cls::allocating_global_matrices_fn();
  main_ns::Matrices_ns::Matrices_Skyline_cls::allocating_local_matrices_fn();
}

/*
###################################################################################################
Purpose: This function allocates global matrices. In this module we consider full matrices.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/18/2018 - Subroutine initiated.


###################################################################################################
*/

void main_ns::Matrices_ns::Matrices_cls::allocating_global_matrices_fn()
{

  std::cout << " -allocating global matrices ..." << std::endl;

  JD = new int[NEqM];
  NTK = new int[NEqM];

  Skyline(NEqM, NEl, NNode, NDOF, NTK, INod, ID, JD);

  K_S = new double[JD[NEqM - 1]]; // Stiffness Matrix
  C_S = new double[JD[NEqM - 1]]; // Damping matrix
  M_S = new double[JD[NEqM - 1]]; // Mass matrix

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

  std::cout << " -allocating DRM matrices ..." << std::endl;
  K_eb = new double *[Model->NDim * Model->NNLayer]; //
  for (int i = 0; i < (Model->NDim * Model->NNLayer); i++)
  {
    K_eb[i] = new double[Model->NDim * Model->NNBndry];
  }

  C_eb = new double *[Model->NDim * Model->NNLayer]; //
  for (int i = 0; i < (Model->NDim * Model->NNLayer); i++)
  {
    C_eb[i] = new double[Model->NDim * Model->NNBndry];
  }

  M_eb = new double *[Model->NDim * Model->NNLayer]; //
  for (int i = 0; i < (Model->NDim * Model->NNLayer); i++)
  {
    M_eb[i] = new double[Model->NDim * Model->NNBndry];
  }

  ND_b = new int[Model->NNBndry * Model->NDim];
  ND_e = new int[Model->NNLayer * Model->NDim];

  // Filling the index for layered nodes
  for (int i = 0; i < Model->NNLayer; i++)
  {
    for (int j = 0; j < Model->NDim; j++)
    {
      ND_e[j * Model->NNLayer + i] = DiscretizedModel->ID[DiscretizedModel->NoLayer_DRM[i]][j];
    }
  }

  // Filling the index for boundary nodes
  for (int i = 0; i < Model->NNBndry; i++)
  {
    for (int j = 0; j < Model->NDim; j++)
    {
      ND_b[j * Model->NNBndry + i] = DiscretizedModel->ID[DiscretizedModel->NoBndry_DRM[i]][j];
    }
  }

  std::cout << " Done with allocation, successfully." << std::endl;
}

/*
###################################################################################################
Purpose: This function assembles local matrices into the full matrices.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/18/2018 - Subroutine initiated.

###################################################################################################
*/

void assemble_local_to_global_fn()
{

  int l, n, i, j, ij;

  for (l = 0; l < NEqEl; l++)
  {
    for (n = 0; n < NEqEl; n++)
    {
      i = ND[l];
      j = ND[n];
      if (i > j)
        continue;
      ij = JD[j] + i - j;

      K_S[ij] = K_S[ij] + Ke[l][n];
      C_S[ij] = C_S[ij] + Ce[l][n];
      M_S[ij] = M_S[ij] + Me[l][n];
    }
  }
}
