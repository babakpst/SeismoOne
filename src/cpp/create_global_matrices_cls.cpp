
#include "../include/create_global_matrices_cls.h"

main_ns::Matrices_ns::Matrices_cls::Matrices_cls
                                (main_ns::discretization_ns::discretization_cls* aDiscretization,
                                 main_ns::model_ns::model_cls* aModel):
                                DiscretizedModel(aDiscretization),
                                Model(aModel)
                                {}

/*
###################################################################################################
Purpose: This function computes the local matrices (mass, damping, and stiffness), and trasmits 
the matrices to be assembled in the global matrices.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 05/14/2018 - Subroutine initiated.
V0.01: 06/02/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/

void main_ns::Matrices_ns::Matrices_cls::compute_elemental_matrices_fn
                                             (int iel, double Rho, double E){

double WX;              // weight in Gauss integration scheme - in the x coordinate
double DJ;              // Jacobian
double DJI;             // Jacobian inverse
double DETJ;            // determinant of Jacobian
double FAC;             // temporary factor to hold the coefficient for numerical integration
double c;               // speed of wave

double * DFX;           // array to hold the value of shape functions at point x

double ** Psi_Psi_T;    // holds the phi*phi^T
double ** PsiX_PsiX_T;

// allocating arrays
DFX = new double [Model->NNode];

Psi_Psi_T = new double*[NEqEl];  // node connectivity
for(int i=0;i<NEqEl;i++){
  Psi_Psi_T[i]=new double[NEqEl];
}

PsiX_PsiX_T = new double*[NEqEl];  // node connectivity
for(int i=0;i<NEqEl;i++){
  PsiX_PsiX_T[i]=new double[NEqEl];
}

// Integrating over the element
  for (int lx=0; lx<Model->NInt; lx++){

    SF->x1  = SF->XInt[lx];
    WX      = SF->WInt[lx];
    
    // Shape functions and differential of shape functions at this local point
    SF->ShapeFunctions();
    SF->DifferentialOfShapeFunctions();

    // Jacobian
    DJ   = 0.0;
      for (int i=0;i<Model->NInt;i++) {
          DJ   += XT[0][i] * SF->DFXI[i];
      }

    DETJ = DJ;   // Jacobian
    FAC  = WX * DETJ;

      if (DETJ <= 0.0) std::cout << "Jacobian is negative!!!" << std::endl;

    // CALCULATING THE INVERSE OF THE JACOBIAN
    DJI = 1.0 / DETJ;

      for (int i=0;i<Model->NInt;i++) {
        DFX [i] = SF->DFXI[i] * DJI;
      }

      for (int i=0;i<Model->NInt;i++) {
        for (int j=0;j<Model->NInt;j++) {
          Psi_Psi_T   [i][j]  = SF->Fn[i] * SF->Fn[j] * FAC;          
          PsiX_PsiX_T [i][j]  = DFX[i] * DFX[j] * FAC;
        }
      }  

      for (int i=0;i<NEqEl;i++) {
        for (int j=0;j<NEqEl;j++) {
          // mass matrix
          Me [i][j] +=  Rho * Psi_Psi_T [i][j] ;

          // stiffness matrix
          Ke [i][j] +=  E   * PsiX_PsiX_T [i][j] ;

          // damping matrix
          Ce [i][j] = 0.0;
        }
      }

  }


c = sqrt(E/Rho);  // speed of wave
if (iel == (Model->NEl)-1) {
 Ce[NEqEl-1][NEqEl-1] = E/c;
}

  for(int i=0;i<NEqEl;i++){
   delete []Psi_Psi_T[i];
  }
delete []Psi_Psi_T;

  for(int i=0;i<NEqEl;i++){
   delete []PsiX_PsiX_T[i];
  }
delete []PsiX_PsiX_T;

}