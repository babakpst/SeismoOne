

#include "../include/Fre_Full.h"
using namespace std;



//***************************************************************************************************************************************************
// Solver: LDLT
//***************************************************************************************************************************************************
void LDLT_Freq ( int& NEqM, double **& K_Eff)
{

// = Local Variables ================================================================================================================================
//int i, j, k, l ;         // Loop index
//int tempI;
double * L;

// ==================== Code ========================================================================================================================
//std::cout << " Factorization ... "<< endl;

L   = new double [ 2*NEqM ];  // Identifications

  for (int j=0; j<(2*NEqM); j++){
    //cout << j << " out of "<< NEqM<< endl;

    //tempI = 5+j;
    //if (tempI> NEqM) tempI = NEqM;

    //for (i=j+1; i<tempI; i++){
    for (int i=j+1; i<(2*NEqM); i++){
      L[i] =  K_Eff[i][j] / K_Eff[j][j];
    }
    for (int k=j+1; k<(2*NEqM); k++){
      for (int l=j+1; l<(2*NEqM); l++){
        K_Eff[k][l] =  K_Eff[k][l] -L[k] * K_Eff[j][l];
      }
    }
    for (int i=j+1; i<(2*NEqM); i++){
      K_Eff[i][j] = L[i] ;
    }
  }
}


//***************************************************************************************************************************************************
// Solver: Substitution
//***************************************************************************************************************************************************
void Substitute_Freq ( int& NEqM, double *& RHS, double **& K_Eff)
{

// = Local Variables ================================================================================================================================

int k, l ;         // Loop index
double temp;
double * L;

// ==================== Code ========================================================================================================================

L   = new double [ 2*NEqM ];  // Identifications

//cout << "Forward" << endl;
  for (int i=0; i< (2*NEqM); i++) {
    temp = 0.0;
      for (int j = 0; j<i; j++){
        temp += K_Eff[i][j] * RHS[j];
      }
    RHS[i]  = RHS[i]-temp ;
  }

  for (int i=0; i<(2*NEqM); i++){
    RHS[i] = RHS[i]/K_Eff[i][i];
  }

//cout << "Backward" << endl;
  for (int i=0; i<(2*NEqM); i++) {

    k = (2*NEqM)-i-1;
    temp = 0.0;
      for (int j=0; j<i; j++) {
       l=(2*NEqM)-j-1;
        temp += K_Eff[l][k] * L[l];
      }
    L[k] = (RHS[k] - temp) ;
  }

  for (int i=0; i<(2*NEqM); i++) {
    RHS[i] = L[i] ;
  }
}


//***************************************************************************************************************************************************
// Calculating the transfer functions 
//***************************************************************************************************************************************************
void Transfer_Full ( double & alpha1, double & alpha2, int & Wave_Type, int& NEqM, double **& M, double **& C, double **& K, double **& PMat, double **& XYZ, int *&ND_e, int *&ND_b, ofstream& TransferFunc )
{

// = Local Variables ================================================================================================================================

// - Integer variables ------------------------------------------------------------------------------------------------------------------------------
//int    i, j ;               // Loop indices
int    NStep;                   // total number of frequency sweep
int    iStep;                   // loop index

// - Real variables ---------------------------------------------------------------------------------------------------------------------------------
double df ;                     // cyclic frequency increment
double minf ;                   // min cyclic frequency
double maxf ;                   // max cyclic frequency

double E;                       // Elastic Modulus
double Rho;                     // density
double c ;                      // wave velocity in the domain where the DRM boundary is
double freq ;                   // cyclicfrequency increments
double omega ;                  // circular frequency increments

double x ;                      // Coordinate of the node
double u_R ;                    // Initialize the values - Not really necessary
double u_I ;                    // Initialize the values - Not really necessary

double Result_R ;               // Transfer functions on the sruface - real part
double Result_I ;               // Transfer functions on the sruface - imaginary part
double Result   ;               // Transfer functions on the sruface - Total motion

const double pi = 3.14159265358979323846;

// - Strings ----------------------------------------------------------------------------------------------------------------------------------------

// - bool -------------------------------------------------------------------------------------------------------------------------------------------

// - integer arrays ---------------------------------------------------------------------------------------------------------------------------------

// - real arrays ------------------------------------------------------------------------------------------------------------------------------------
double * U;                     // Solution in the frequency domain - The first half of the vector is the real part and the second half is the imaginary part
double * RHS;                   // Right Hand Side of the equation (load vector)  - The first half of the vector is the real part and the second half is the imaginary part

double * U_b_R;                 // The real part of the analytical solution on the boundary
double * U_b_I;                 // The imaginary part of the analytical solution on the boundary
double * U_e_R;                 // The real part of the analytical solution on the layer
double * U_e_I;                 // The imaginary part of the analytical solution on the boundary

double * F_b_R;                 // The real part of the DRM load on the DRM boundary
double * F_b_I;                 // The imaginary part of the DRM load on the DRM boundary
double * F_e_R;                 // The real part of the DRM load on the DRM layer
double * F_e_I;                 // The imaginary part of the DRM load on the DRM layer

double ** K_Eff;                // complex Effective stiffness (real + imaginary) - See notes
double ** K_Eff_eb;             // The effective striffness DRM matrix
double ** C_Eff_eb;             // The effective damping DRM matrix

// - data structure ---------------------------------------------------------------------------------------------------------------------------------
// ==================== Code ========================================================================================================================

cout << "<<<<   Transfer function void   >>>>" << endl;

// Define arrays
cout << "Create arrays ..." << endl;

U     = new double [2*NEqM];     // Solution for each frequency (real part + imaginary part)
RHS   = new double [2*NEqM];     // Solution for each frequency (real part + imaginary part)

K_Eff = new double *[ 2*NEqM ];  // Effective stiffness
  for(int i=0;i<(2*NEqM);i++){
    K_Eff[i]=new double[2*NEqM];
  }

K_Eff_eb  = new double *[ NDim * NNLayer ];  // The effective striffness DRM matrix
  for(int i=0;i<(NDim * NNLayer);i++){
    K_Eff_eb[i]=new double[ NDim * NNBndry ];
  }

C_Eff_eb  = new double *[ NDim * NNLayer ];  // The effective striffness DRM matrix
  for(int i=0;i<(NDim * NNLayer);i++){
    C_Eff_eb[i]=new double[ NDim * NNBndry ];
  }


// Defining the required vectors
U_b_R = new double [ NNBndry * NDim ];
U_b_I = new double [ NNBndry * NDim ];

F_e_R   = new double [ NNLayer * NDim ];
F_e_I   = new double [ NNLayer * NDim ];

// Defining the required vectors
U_e_R = new double [ NNLayer * NDim ];
U_e_I = new double [ NNLayer * NDim ];

F_b_R = new double [ NNBndry * NDim ];
F_b_I = new double [ NNBndry * NDim ];


// Define frequencies
minf = 0.0; 
maxf = alpha1;
df   = alpha2;
 
// material properties of the domain where the DRM boundary is 
E   = PMat[0][0];
Rho = PMat[0][1];
c = sqrt (E/Rho);

// NSteps
NStep = (int)((maxf-minf)/df);  // Total number of steps to cover the frequency range

  for (iStep =1; iStep<=NStep; iStep++){

    freq  = minf + iStep * df;
    omega = 2 * pi * freq ;
    std::cout << "Step: " << iStep << " Frequency: " << freq << endl;

    // Effective stiffness matrix for this frequency
    //std::cout << "Effective stiffness ..." << endl;
      for (int i=0; i<NEqM; i++) {
        for (int j=0; j<NEqM; j++) {
          K_Eff[i][j]           = - omega* omega * M[i][j] + K[i][j];

          K_Eff[i+NEqM][j]      = - omega * C[i][j] ;

          K_Eff[i][j+NEqM]      = - omega * C[i][j] ;

          K_Eff[i+NEqM][j+NEqM] = + omega* omega * M[i][j] - K[i][j];
        }
      }

      // Coefficient matrices for DRM loads for this frequency
      for (int i=0; i<(NDim * NNLayer); i++) {
        for (int j=0;j<(NDim * NNBndry);j++) {
          K_Eff_eb[i][j] = + omega * omega * M_eb[i][j] - K_eb[i][j] ; 
          C_Eff_eb[i][j] = + omega * C_eb[i][j] ; 
        }
      }

      // Initializing displacement and the load vector (RHS)
      for (int i=0; i<(2*NEqM); i++) {
        U[i]    = 0.0 ;
        RHS[i]  = 0.0 ;
      }

    // ------- Working on F_e -----------------------------
      for (int i=0; i<NNBndry * NDim; i++ ) {
        U_b_R [i] = 0.0;
        U_b_I [i] = 0.0;
      }

      for (int i=0; i<NNLayer * NDim; i++ ) {
        F_e_R [i] = 0.0;
        F_e_I [i] = 0.0;
      }

      // Computing the analytical solution at the DRM boundary and the DRM layer nodes (Flat homogeneous domain)
      for (int i=0;i<NNBndry;i++) {
        for (int j=0;j<NDim;j++) {
          x  = XYZ [ NoBndry_DRM[i] ][j];  // Coordinate of the node
          u_R =  0.0;  // Initialize the values - Not really necessary
          u_I =  0.0;  // Initialize the values - Not really necessary
            // Computing the analytical solution - Comment: the one-dimensional wave-motion is identical for both SV and P waves.
            if      ( Wave_Type == 0 ) DRM_PointValues_Freq ( amplitude, x, c, omega, u_R, u_I ); // SV wave
            else if ( Wave_Type == 1 ) DRM_PointValues_Freq ( amplitude, x, c, omega, u_R, u_I ); // P wave - Basically, the same as the SV wave in the one-dimensional simulation

          // Filling the analytical solution vector
          U_b_R  [i * NNBndry * NDim + j] = u_R ;
          U_b_I  [i * NNBndry * NDim + j] = u_I ;
        }
      }

      // 
      for (int i=0;i<NNLayer*NDim;i++) {
        for (int j=0;j<NNBndry*NDim;j++) {
          F_e_R [i] += -K_Eff_eb[i][j] * U_b_R [j] - C_Eff_eb[i][j] * U_b_I [j];
          F_e_I [i] += -K_Eff_eb[i][j] * U_b_I [j] + C_Eff_eb[i][j] * U_b_R [j]; 
        }
      }

    // Assemble the load vector 
      for (int i=0; i <NNLayer*NDim; i++){
        RHS[ ND_e [i]        ] += +F_e_R [i];
        RHS[ ND_e [i] + NEqM ] += -F_e_I [i];
      }
    
    // ------- Working on F_b -----------------------------
      for (int i=0; i<NNLayer * NDim; i++ ) {
        U_e_R [i] = 0.0;
        U_e_I [i] = 0.0;
      }

      for (int i=0; i<NNBndry * NDim; i++ ) {
        F_b_R [i] = 0.0;
        F_b_I [i] = 0.0;
      }

      // Loop on the nodes on the DRM layer to find out the analytical solution (In this case only two nodes)
      for (int i=0;i<NNLayer; i++) {
        for (int j=0;j<NDim;j++) {
          x  = XYZ [ NoLayer_DRM[i] ][j];  // Coordinate of the node
          u_R =  0.0;  // Initialize the values - Not really necessary
          u_I =  0.0;  // Initialize the values - Not really necessary
            // Computing the analytical solution - Comment: the one-dimensional wave-motion is identical for both SV and P waves.
            if      ( Wave_Type == 0 ) DRM_PointValues_Freq ( amplitude, x, c, omega, u_R, u_I ); // SV wave
            else if ( Wave_Type == 1 ) DRM_PointValues_Freq ( amplitude, x, c, omega, u_R, u_I ); // P wave - Basically, the same as the SV wave in the one-dimensional simulation

          // Filling the analytical solution vector
          U_e_R  [i * NNBndry * NDim + j] = u_R ;
          U_e_I  [i * NNBndry * NDim + j] = u_I ;
        }
      }

      // Multiply the Mass, Damp, and Stiffness matrix by the vector
      for (int i=0;i<NNBndry*NDim;i++) {
        for (int j=0;j<NNLayer*NDim;j++) {
          F_b_R [i] += +K_Eff_eb[j][i] * U_e_R [j] + C_Eff_eb[j][i] * U_e_I [j];
          F_b_I [i] += +K_Eff_eb[j][i] * U_e_I [j] - C_Eff_eb[i][j] * U_e_R [j]; 
        }
      }

      // Assemble the load vector 
      for (int i=0; i <NNBndry*NDim; i++){
        RHS [ ND_b [i]        ] += +F_b_R[i]; 
        RHS [ ND_b [i] + NEqM ] += -F_b_I[i]; 
      }

    // Reduction the coefficient matrix ( Effective Stiffness Matrix )
    //Reduce_Full (NEqM, K, Check);
    LDLT_Freq ( NEqM,  K_Eff);

    // SOLVE
    //Gauss_El_Full ( NEqM, UN, K);
    //Gaussian ( NEqM, UN, K);
    Substitute_Freq ( NEqM, RHS, K_Eff);

    // time history of the solution at some particular nodes
    Result_R = RHS[ NEqM - 1];
    Result_I = RHS[ NEqM * 2 - 1 ] ;
    Result   = sqrt( Result_R*Result_R + Result_I*Result_I );

    TransferFunc << setw(20) << freq << setw(20) << Result_R  << setw(20) << Result_I  << setw(20) << Result << endl;

  }

// Deallocate vectors
delete U_b_R ;
delete U_b_I ;
delete F_e_R ;
delete F_e_I ;

// Deallocate vectors
delete U_e_R ;
delete U_e_I ;
delete F_b_R ;
delete F_b_I ;

}
