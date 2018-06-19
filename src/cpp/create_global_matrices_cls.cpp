
#include "../include/create_global_matrices_cls.h"

main_ns::Matrices_ns::Matrices_cls::Matrices_cls
                                (main_ns::discretization_ns::discretization_cls* aDiscretization,
                                 main_ns::model_ns::model_cls* aModel):
                                DiscretizedModel(aDiscretization),
                                Model(aModel)
                                {}

/*
###################################################################################################
Purpose: This function allocates local matrices for each element.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 05/14/2018 - Subroutine initiated.
V0.01: 05/15/2018 - Initiated: Compiled without error for the first time.
V1.00: 06/18/2018 - 

###################################################################################################
*/

void main_ns::Matrices_ns::Matrices_Full_cls::allocating_local_matrices_fn(){

NEqEl  = Model->NDOF * Model->NNode ;

// Allocating Elemental Matrices
Ke = new double*[NEqEl];
for(int i=0;i<NEqEl;i++){
  Ke[i]=new double[NEqEl];
}

Ce = new double*[NEqEl];
for(int i=0;i<NEqEl;i++){
  Ce[i]=new double[NEqEl];
}

Me = new double*[NEqEl];
for(int i=0;i<NEqEl;i++){
  Me[i]=new double[NEqEl];
}

XT = new double*[Model->NDim];
for(int i=0;i<Model->NNode;i++){
  XT[i]=new double[Model->NNode];
}

Fe  = new double [NEqEl];
ND = new int [NEqEl];

}


/*
###################################################################################################
Purpose: This function computes the local matrices for each element and assembles the local 
matrices into the golabal matrices.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 05/14/2018 - Subroutine initiated.
V0.01: 05/15/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/

void main_ns::Matrices_ns::Matrices_cls::assembling_local_matrices_into_global_matrices_fn(){

    int ElementPercent;     // To show the progress in assembly
    double AssemblyPercentage; //

// defining the shape function, SF, depending on the order of the shape function
if (Model->OrderOfShapeFunc == 1){
   SF = new main_ns::ShapeFunctions_ns::ShapeFunctions_FirstOrder_cls(Model->NInt, Model->NNode);}
else if (Model->OrderOfShapeFunc == 2){
   SF = new main_ns::ShapeFunctions_ns::ShapeFunctions_SecondOrder_cls(Model->NInt, Model->NNode);}

SF->Retrieving_Gauss_Points_fn(); // Extracting quadratures

// In order to print the progress in computing local matrices, in 10 steps,  we define this var.
ElementPercent = (int)( Model->NEl/10.0);

  // computing element matrices and assembling them
  for (int iel=0; iel<Model->NEl; iel++){

    // writing down the percentage of the work done on screen
    if ((iel % ElementPercent) == 0) {
      AssemblyPercentage = ((double) iel/Model->NEl)*100.0;
      std::cout << "Assembly progress:  %" << AssemblyPercentage << std::endl;
    }
    
      // extracting the coordinates of this element
      for (int i=0; i<Model->NDim; i++) {
        for (int j=0; j<Model->NNode; j++) {
          XT [i][j] =  DiscretizedModel->XYZ [ DiscretizedModel->INod[j][iel] ][i]; 
        }
      }

    // Initializing element matrices (stiffness, damping, and mass), and the load vector
      for (int i=0; i<NEqEl; i++) {
          for (int j=0; j<NEqEl; j++) {
            Ke[i][j] = 0.0;
            Ce[i][j] = 0.0;
            Me[i][j] = 0.0;
          }
        Fe[i] = 0.0;
      }

    // Material Property of this element
    MType = DiscretizedModel->MTel [iel];    // Material property type
    E   = Model->PMat[MType][0];  // Elastic modulus of this material
    Rho = Model->PMat[MType][1];  // Density of this material 

    
    compute_elemental_matrices_fn(iel, Rho, E);

      for (int i=0;i<Model->NNode;i++) {
        for (int j=0;j<Model->NDOF;j++) {
          ND [ j * Model->NNode + i ] = DiscretizedModel->ID [DiscretizedModel->INod [i][iel]][j];
        }
      }

    // assemble local mass, stiffness, and damping matrices in the global matrices
    assemble_local_to_global_fn();
  }

// Deallocating Element Matrices
for(int i=0;i<NEqEl;i++){
 delete []Ke[i];
}
delete []Ke;

for(int i=0;i<NEqEl;i++){
 delete []Ce[i];
}
delete []Ce;

for(int i=0;i<NEqEl;i++){
 delete []Me[i];
}
delete []Me;

//
//for(int i=0;i<NDim;i++){
// delete []XT[i];
//}
//delete []XT;
//

delete Fe;
delete ND;

}

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