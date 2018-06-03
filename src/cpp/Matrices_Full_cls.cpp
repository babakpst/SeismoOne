

#include "../include/Matrices_Full_cls.h"


// Constructor
main_ns::Matrices_Full_ns::Matrices_Full_cls::Matrices_Full_cls
                            (main_ns::discretization_ns::discretization_cls* aDiscretization,
                             main_ns::model_ns::model_cls* aModel):
                             main_ns::Matrices_ns::Matrices_cls(aDiscretization,aModel){
                                  allocating_global_matrices_fn(); // allocating global matrices
                                  allocating_local_matrices_fn(); // allocating local matrices matrices
                             }


// 
/*
###################################################################################################
Purpose: This function Allocating global matrices. In this module we consider full matrices.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 05/14/2018 - Subroutine initiated.
V0.01: 05/15/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/

void main_ns::Matrices_Full_ns::Matrices_Full_cls::allocating_global_matrices_fn(){

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

// allocating local matrices for each element
void main_ns::Matrices_Full_ns::Matrices_Full_cls::allocating_local_matrices_fn(){

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

void compute_elemental_matrices_fn();

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

void main_ns::Matrices_Full_ns::Matrices_Full_cls::
     assembling_local_matrices_into_global_matrices_fn(){

main_ns::ShapeFunctions_ns::ShapeFunctions_cls ShapeFunctions;
ShapeFunctions.Retrieving_Gauss_Points_fn(); // Extracting quadratures

// In order to print the progress in computing local matrices, in 10 steps,  we define this var.
int ElementPercent = (int)( Model->NEl/10.0);

  // computing element matrices and assembling them
  for (unsigned int iel=0; iel<Model->NEl; iel++) {
    // This if condition prints the work done on screen
    if (( iel % ElementPercent ) == 0) {
      AssemblyPercentage = ((double) iel/NEl)*(double)100;
      std::cout << "Assembly progress:  %" << AssemblyPercentage << sted::endl;
    }

    // Material Property of this element
    MType = MTel [iel];    // Material property type
    E   = PMat[MType][0];  // Elastic modulus of this material
    Rho = PMat[MType][1];  // Density of this material 

      // coordinates of this element
      for (unsigned int i=0; i<Model->NDim; i++) {
        for (unsigned int j=0; j<Model->NNode; j++) {
          XT [i][j] =  XYZ [ INod[j][iel] ][i]; 
        }
      }

    // Initializing element matrices (stiffness, damping, and mass), and the load vector
      for (unsigned int i=0; i<Model->NEqEl; i++) {
          for (unsigned int j=0; j<Model->NEqEl; j++) {
            Ke[i][j] = 0.0;
            Ce[i][j] = 0.0;
            Me[i][j] = 0.0;
          }
        Fe[i] = 0.0;
      }



    // computing elemental matrices
    if (OShFunc == 1)
      MassDampStiffS_1D_first_Full ( iel, NEl, NInt, NNode, NEqEl, Rho, E, XT, Me, Ce, Ke, Gauss_PNT);
    else if (OShFunc == 2)
      MassDampStiffS_1D_second_Full ( iel, NInt, NNode, NEqEl, Rho, E, XT, Me, Ce, Ke, Gauss_PNT);

    
    // force vector of the element
    // not applicable in this code

    // assemble mass, stiffness, damping and force matrices
      for (int i=0;i<NNode;i++) {
        for (int j=0;j<NDOF;j++) {
          ND [ j * NNode + i ] = ID [INod [i][iel]][j];
        }
      }


/*
    // Check element matrices
    Check << "Damping matrix" << "\n";
      for (int i = 0; i<NEqEl; i++){
        for (int j = 0; j<NEqEl; j++){
          Check << setw(20) << Ce[i][j] ;
        }
        Check <<  endl;
      }

    Check <<  endl;
    Check <<  endl;
    Check <<  endl;
*/

    AssembleMassDampingStiffForceFull (NEqEl, ND, Ke, Ce, Me, K, C, M);



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

/*
for(int i=0;i<NDim;i++){
 delete []XT[i];
}
delete []XT;
*/

delete Fe;
delete ND;



}



//***************************************************************************************************************************************************
// computing element matrices - first order
//***************************************************************************************************************************************************
void MassDampStiffS_1D_first_Full ( int& iel, int& NEl, int& NInt, int& NNode, int& NEqEl, double& Rho, double& E, double **& XT, double ** &Me, double ** &Ce,  double ** &Ke, Gauss& Gauss_PNT)
{

// = Local Variables ================================================================================================================================

// - Integer variables ------------------------------------------------------------------------------------------------------------------------------

// - Real variables ---------------------------------------------------------------------------------------------------------------------------------
double WX;              // weight in Gauss integration scheme
double DJ;              // Jacobian
double DJI;             // Jacobian inverse
double DETJ;            // determinant of Jacobian
double FAC;             // temporary factor
double c;               // speed of wave

// - Strings ----------------------------------------------------------------------------------------------------------------------------------------

// - bool -------------------------------------------------------------------------------------------------------------------------------------------

// - real arrays ------------------------------------------------------------------------------------------------------------------------------------
double * DFX;

double ** Psi_Psi_T; 
double ** PsiX_PsiX_T ;

// - data structure ---------------------------------------------------------------------------------------------------------------------------------
SF_1D_2N    SF;
DSF_1D_2N  DSF;

// ==================== Code ========================================================================================================================

// allocating arrays
DFX = new double [NNode];

Psi_Psi_T = new double*[NEqEl];  // node connectivity
for(int i=0;i<NEqEl;i++){
  Psi_Psi_T[i]=new double[NEqEl];
}

PsiX_PsiX_T = new double*[NEqEl];  // node connectivity
for(int i=0;i<NEqEl;i++){
  PsiX_PsiX_T[i]=new double[NEqEl];
}

// Integrating over the element
  for (int lx=0;lx<NInt;lx++) {

    SF.X1  = Gauss_PNT.XInt [ lx ] ;
    DSF.X1 = Gauss_PNT.XInt [ lx ] ;
    WX     = Gauss_PNT.WInt [ lx ] ;

    // SHAPE FUNCTIONS AND DIFFERENTIAL OF SHAPE FUNCTIONS
    ShapeFunc_1D_2N (  SF.X1, SF.FN[0], SF.FN[1] ) ;
    Dif_ShapeFunc_1D_2N ( DSF.DFXI[0], DSF.DFXI[1] ) ;


    // Jacobian
    DJ   = XT[0][0] * DSF.DFXI[0] + XT[0][1] * DSF.DFXI[1];


    DETJ = DJ ;   // Jacobian
    FAC  = WX * DETJ ;

      if ( DETJ <= 0.0 ) std::cout << "Jacobian is negative!!!" << endl;

    // CALCULATING THE INVERSE OF THE JACOBIAN
    DJI = 1.0 / DETJ ;

    DFX [0] = DSF.DFXI[0] * DJI ;
    DFX [1] = DSF.DFXI[1] * DJI ;

    Psi_Psi_T [0][0]  = SF.FN[0] * SF.FN[0] * FAC ;
    Psi_Psi_T [0][1]  = SF.FN[0] * SF.FN[1] * FAC ;
    Psi_Psi_T [1][0]  = SF.FN[1] * SF.FN[0] * FAC ;
    Psi_Psi_T [1][1]  = SF.FN[1] * SF.FN[1] * FAC ;

    PsiX_PsiX_T [0][0]  = DFX[0] * DFX[0] * FAC ;
    PsiX_PsiX_T [0][1]  = DFX[0] * DFX[1] * FAC ;
    PsiX_PsiX_T [1][0]  = DFX[1] * DFX[0] * FAC ;
    PsiX_PsiX_T [1][1]  = DFX[1] * DFX[1] * FAC ;

      for (int i=0;i<NEqEl;i++) {
        for (int j=0;j<NEqEl;j++) {
          // mass matrix
          Me [i][j] = Me [i][j] + Rho * Psi_Psi_T [i][j] ;

          // stiffness matrix
          Ke [i][j] = Ke [i][j] + E   * PsiX_PsiX_T [i][j] ;

          // damping matrix
          Ce [i][j] = 0.0;

        }
      }

  }


c = sqrt(E/Rho);  // speed of wave
if (iel == NEl-1) {
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


//***************************************************************************************************************************************************
// computing element matrices - second order
//***************************************************************************************************************************************************
void MassDampStiffS_1D_second_Full ( int& iel, int& NInt, int& NNode, int& NEqEl, double& Rho, double& E, double **& XT, double ** &Me, double ** &Ce,  double ** &Ke, Gauss& Gauss_PNT)
{

// = Local Variables ================================================================================================================================

// - Integer variables ------------------------------------------------------------------------------------------------------------------------------

// - Real variables ---------------------------------------------------------------------------------------------------------------------------------
double WX;              // weight in Gauss integration scheme
double DJ;              // Jacobian
double DJI;             // Jacobian inverse
double DETJ;            // determinant of Jacobian
double FAC;             // temporary factor
double c;               // speed of wave

// - Strings ----------------------------------------------------------------------------------------------------------------------------------------

// - bool -------------------------------------------------------------------------------------------------------------------------------------------

// - real arrays ------------------------------------------------------------------------------------------------------------------------------------
double * DFX;

double ** Psi_Psi_T; 
double ** PsiX_PsiX_T ;

// - data structure ---------------------------------------------------------------------------------------------------------------------------------
SF_1D_3N    SF;
DSF_1D_3N  DSF;

// ==================== Code ========================================================================================================================

// allocating arrays
DFX = new double [NNode];

Psi_Psi_T = new double*[NEqEl];  // node connectivity
for(int i=0;i<NEqEl;i++){
  Psi_Psi_T[i]=new double[NEqEl];
}

PsiX_PsiX_T = new double*[NEqEl];  // node connectivity
for(int i=0;i<NEqEl;i++){
  PsiX_PsiX_T[i]=new double[NEqEl];
}

// Integrating over the element
  for (int lx=0;lx<NInt;lx++) {

    SF.X1  = Gauss_PNT.XInt [ lx ] ;
    DSF.X1 = Gauss_PNT.XInt [ lx ] ;
    WX     = Gauss_PNT.WInt [ lx ] ;

    // SHAPE FUNCTIONS AND DIFFERENTIAL OF SHAPE FUNCTIONS
    ShapeFunc_1D_3N (  SF.X1, SF.FN[0], SF.FN[1], SF.FN[2] ) ;
    Dif_ShapeFunc_1D_3N ( DSF.X1, DSF.DFXI[0], DSF.DFXI[1], DSF.DFXI[2] ) ;

    // Jacobian
    DJ   = XT[0][0] * DSF.DFXI[0] + XT[0][1] * DSF.DFXI[1] + XT[0][2] * DSF.DFXI[2];

    DETJ = DJ ;   // Jacobian
    FAC  = WX * DETJ ;

      if ( DETJ <= 0.0 ) std::cout << "Jacobian is negative!!!" << endl;

    // CALCULATING THE INVERSE OF THE JACOBIAN
    DJI = 1.0 / DETJ ;

    DFX [0] = DSF.DFXI[0] * DJI ;
    DFX [1] = DSF.DFXI[1] * DJI ;
    DFX [2] = DSF.DFXI[2] * DJI ;


    Psi_Psi_T [0][0]  = SF.FN[0] * SF.FN[0] * FAC ;
    Psi_Psi_T [0][1]  = SF.FN[0] * SF.FN[1] * FAC ;
    Psi_Psi_T [0][2]  = SF.FN[0] * SF.FN[2] * FAC ;

    Psi_Psi_T [1][0]  = SF.FN[1] * SF.FN[0] * FAC ;
    Psi_Psi_T [1][1]  = SF.FN[1] * SF.FN[1] * FAC ;
    Psi_Psi_T [1][2]  = SF.FN[1] * SF.FN[2] * FAC ;

    Psi_Psi_T [2][0]  = SF.FN[2] * SF.FN[0] * FAC ;
    Psi_Psi_T [2][1]  = SF.FN[2] * SF.FN[1] * FAC ;
    Psi_Psi_T [2][2]  = SF.FN[2] * SF.FN[2] * FAC ;

    PsiX_PsiX_T [0][0]  = DFX[0] * DFX[0] * FAC ;
    PsiX_PsiX_T [0][1]  = DFX[0] * DFX[1] * FAC ;
    PsiX_PsiX_T [0][2]  = DFX[0] * DFX[2] * FAC ;

    PsiX_PsiX_T [1][0]  = DFX[1] * DFX[0] * FAC ;
    PsiX_PsiX_T [1][1]  = DFX[1] * DFX[1] * FAC ;
    PsiX_PsiX_T [1][2]  = DFX[1] * DFX[2] * FAC ;

    PsiX_PsiX_T [2][0]  = DFX[2] * DFX[0] * FAC ;
    PsiX_PsiX_T [2][1]  = DFX[2] * DFX[1] * FAC ;
    PsiX_PsiX_T [2][2]  = DFX[2] * DFX[2] * FAC ;


      for (int i=0;i<NEqEl;i++) {
        for (int j=0;j<NEqEl;j++) {
          // mass matrix
          Me [i][j] = Me [i][j] + Rho * Psi_Psi_T [i][j] ;

          // stiffness matrix
          Ke [i][j] = Ke [i][j] + E   * PsiX_PsiX_T [i][j] ;

          // damping matrix
          Ce [i][j] = 0.0;

        }
      }

  }


c = sqrt(E/Rho);  // speed of wave
if (iel == 0) {
 Ce[0][0] = E/c;
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







