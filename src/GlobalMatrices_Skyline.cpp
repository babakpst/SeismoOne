
#include "GlobalMatrices_Skyline.h"
using namespace std;


//***************************************************************************************************************************************************
// computing global matrices
//***************************************************************************************************************************************************
void GlobalMatrices_Skyline ( int & NInt, int& NDim, int& NDOF, int& NNode, int *& MTel, int **& INod, int **& ID, double **& XYZ, double *& M_S, double *& C_S, double *& K_S, int *& JD) 
{

// = Local Variables ================================================================================================================================

// - Integer variables ------------------------------------------------------------------------------------------------------------------------------
int MType;              // Material type
int NEqEl;              // Number of equations of each element
int ElementPercent;     // To show the progress in assembly
// - Real variables ---------------------------------------------------------------------------------------------------------------------------------
double E;               // elastic modulus
double Rho;             // density
double AssemblyPercentage; //

// - Strings ----------------------------------------------------------------------------------------------------------------------------------------

// - bool -------------------------------------------------------------------------------------------------------------------------------------------

// - integer arrays ---------------------------------------------------------------------------------------------------------------------------------
int    * ND;            // element constraints

// - real arrays ------------------------------------------------------------------------------------------------------------------------------------
double * Fe;            // element force vector

double ** XT;           // Array to store coordinates of the element
double ** Ke;           // element stiffness matrix
double ** Ce;           // element damping matrix
double ** Me;           // element mass matrix

// - data structure ---------------------------------------------------------------------------------------------------------------------------------
Gauss   Gauss_PNT ; // Defining Gauss points

// ==================== Code ========================================================================================================================

cout << "Global Matrices - Skyline method" << endl; 

NEqEl  = NDOF * NNode ;
GAUSS_Quad_POINTS ( NInt, Gauss_PNT ) ;
ElementPercent = (int)(NEl / 10.0);

// Allocating Elemental Matrices

Me = new double*[NEqEl];
for(int i=0;i<NEqEl;i++){
  Me[i]=new double[NEqEl];
}


Ce = new double*[NEqEl];
for(int i=0;i<NEqEl;i++){
  Ce[i]=new double[NEqEl];
}


Ke = new double*[NEqEl];
for(int i=0;i<NEqEl;i++){
  Ke[i]=new double[NEqEl];
}


XT = new double*[NDim];
for(int i=0;i<NNode;i++){
  XT[i]=new double[NNode];
}



Fe  = new double [NEqEl];
ND  = new int [NEqEl];

  // computing element matrices and assembling them
  for (int iel=0;iel<NEl;iel++) {

    if (( iel % ElementPercent ) == 0) {
      AssemblyPercentage = ((double) iel/NEl)*(double)100;
      std::cout << "Assembly progress:  %" << AssemblyPercentage << endl;
    }

    MType = MTel [iel]; // material property of this element

    // Initialize element matrices
      for (int i=0; i<NEqEl; i++) {
          for (int j=0; j<NEqEl; j++) {
            Ke[i][j] = 0.0;
            Ce[i][j] = 0.0;
            Me[i][j] = 0.0;
          }
        Fe[i] = 0.0;
      }



    // Material Property of this element
    E   = PMat[MType][0];
    Rho = PMat[MType][1];



      // coordinates of the element
      for (int i=0; i<NDim; i++) {
        for (int j=0; j<NNode; j++) {
          XT [i][j] =  XYZ [ INod[j][iel] ][i]; 
        }
      }

    // computing elemental matrices
    if (OShFunc == 1)
      MassDampStiffS_1D_first_Skyline ( iel, NEl, NInt, NNode, NEqEl, Rho, E, XT, Me, Ce, Ke, Gauss_PNT);
    else if (OShFunc == 2)
      MassDampStiffS_1D_second_Skyline ( iel, NInt, NNode, NEqEl, Rho, E, XT, Me, Ce, Ke, Gauss_PNT);

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
    Check << "Mass matrix" << "\n";
      for (int i = 0; i<NEqEl; i++){
        for (int j = 0; j<NEqEl; j++){
          Check << setw(20) << Me[i][j] ;
        }
        Check <<  endl;
      }


    Check <<  endl;
    Check <<  endl;
    Check <<  endl;
*/




//for (int i =0; i<NEqEl; i++){
//for (int j =0; j<NEqEl; j++){
//cout << Me[i][j] << "   " << Ce[i][j] << "   " << Ke[i][j] << endl; 
//}
//}
    AssembleMassDampingStiffForceSkyline ( NEqEl, ND, Ke, Ce, Me, K_S, C_S, M_S, JD );

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



for(int i=0;i<NDim;i++){
  delete []XT[i];
  }
delete []XT;

delete Fe;
delete ND;


std::cout << "====== End Assembly - Skyline method ======" << endl; 

}


//***************************************************************************************************************************************************
// computing element matrices - first order
//***************************************************************************************************************************************************
void MassDampStiffS_1D_first_Skyline ( int& iel, int& NEl, int& NInt, int& NNode, int& NEqEl, double& Rho, double& E, double **& XT, double ** &Me, double ** &Ce,  double ** &Ke, Gauss& Gauss_PNT)
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
void MassDampStiffS_1D_second_Skyline ( int& iel, int& NInt, int& NNode, int& NEqEl, double& Rho, double& E, double **& XT, double ** &Me, double ** &Ce,  double ** &Ke, Gauss& Gauss_PNT)
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
void DRM_Matrices_Skyline( int & NNBndry, int &NNLayer, double *& K_S, double *& C_S, double *& M_S, double **& K_eb, double **& C_eb, double **& M_eb, int *&ND_e, int *&ND_b , int *& JD)  
{

int ij,i,j,l,n;   // Loop indices

// - Code ---------------------------------------------------------------------
std::cout << "Create DRM matrices ..." << endl; 

  for ( l = 0; l < NNLayer * NDim; l++) {
    for ( n = 0; n < NNBndry * NDim; n++) {

      i = ND_e [ l ] ;
      j = ND_b [ n ] ;


      if (i>j) {
        ij = i;
        i  = j;
        j  = ij; 
        std::cout << " i and j replaced"<< endl;
      }

      ij = JD[j]+i-j;

//cout<< "DRM" << ij << " K_S  "<< K_S[ ij ] <<  " M_S "<< M_S[ ij ] << " C_S  "<< C_S[ ij ] <<  endl;
   
      K_eb[l][n] = K_S[ ij ];
      C_eb[l][n] = C_S[ ij ];
      M_eb[l][n] = M_S[ ij ];

//      cin.get();
    }
  }

std::cout << "DRM matrices created." << endl; 


}








