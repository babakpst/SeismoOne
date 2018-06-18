#include "../include/create_full_matrices_cls.h"

// Constructor: we also create and allocate matrices
main_ns::Matrices_ns::Matrices_Full_cls::Matrices_Full_cls
                                  (main_ns::discretization_ns::discretization_cls* aDiscretization,
                                   main_ns::model_ns::model_cls* aModel):
                                   main_ns::Matrices_ns::Matrices_cls(aDiscretization,aModel){
                    main_ns::Matrices_ns::Matrices_Full_cls::allocating_global_matrices_fn();
                    main_ns::Matrices_ns::Matrices_Full_cls::allocating_local_matrices_fn(); 
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

void main_ns::Matrices_ns::Matrices_Full_cls::allocating_global_matrices_fn(){

std::cout<< " -allocating global matrices ..." << std::endl;
K  = new double *[DiscretizedModel->NEqM];  // Stiffness Matrix
  for(int i=0;i<DiscretizedModel->NEqM;i++){
    K[i]=new double[DiscretizedModel->NEqM];
  }

C  = new double *[DiscretizedModel->NEqM];  // Damping matrix
  for(int i=0;i<DiscretizedModel->NEqM;i++){
    C[i]=new double[DiscretizedModel->NEqM];
  }

M  = new double *[DiscretizedModel->NEqM];  // Mass matrix
  for(int i=0;i<DiscretizedModel->NEqM;i++){
    M[i]=new double[DiscretizedModel->NEqM];
  }

F  = new double [DiscretizedModel->NEqM] ;

std::cout << " -initializing global matrices ..." << std::endl;
for (int i=0; i<DiscretizedModel->NEqM; i++) {
    for (int j=0; j<DiscretizedModel->NEqM; j++) {
      M[i][j] = 0.0;
      C[i][j] = 0.0;
      K[i][j] = 0.0;
    }
  F[i]=0.0;
}


std::cout<< " -allocating DRM matrices ..." << std::endl;
K_eb  = new double *[Model->NDim * Model->NNLayer];  // 
  for(int i=0;i<(Model->NDim * Model->NNLayer);i++){
    K_eb[i]=new double[ Model->NDim * Model->NNBndry ];
  }

C_eb  = new double *[ Model->NDim * Model->NNLayer ];  // 
  for(int i=0;i<(Model->NDim * Model->NNLayer);i++){
    C_eb[i]=new double[ Model->NDim * Model->NNBndry ];
  }

M_eb  = new double *[ Model->NDim * Model->NNLayer ];  // 
  for(int i=0; i<(Model->NDim * Model->NNLayer); i++){
    M_eb[i]= new double[ Model->NDim * Model->NNBndry ];
  }


ND_b = new int  [ Model->NNBndry * Model->NDim ];
ND_e = new int  [ Model->NNLayer * Model->NDim ];

// Filling the index for layered nodes
  for (int i=0;i<Model->NNLayer;i++) {
    for ( int j=0;j<Model->NDim;j++) {
      ND_e [ j * Model->NNLayer + i ] = DiscretizedModel->ID [ DiscretizedModel->NoLayer_DRM [ i ] ][j];
    }
  }

  // Filling the index for boundary nodes
  for ( int i=0;i<Model->NNBndry;i++) {
    for (int j=0;j<Model->NDim;j++) {
      ND_b [ j * Model->NNBndry + i ] = DiscretizedModel->ID [ DiscretizedModel->NoBndry_DRM[i]][j];
    }

  }

std::cout<< " Done with allocation, successfully." << std::endl;
}


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

void main_ns::Matrices_ns::Matrices_Full_cls::assembling_local_matrices_into_global_matrices_fn(){

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
Purpose: This function assembles local matrices into the full matrices.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/18/2018 - Subroutine initiated.

###################################################################################################
*/

void assemble_local_to_global_fn(){

int L;
int N;

  for (int i=0; i<NEqEl; i++) {
    for (int j=0; j<NEqEl; j++) {
      L = ND [i] ;
      N = ND [j] ;
        if (L == -1 || N == -1) continue;
      K [L][N] = K [L][N] + Ke [i][j];
      C [L][N] = C [L][N] + Ce [i][j];
      M [L][N] = M [L][N] + Me [i][j];
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





