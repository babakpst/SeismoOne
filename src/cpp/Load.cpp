

#include "../include/Load.h"

//***************************************************************************************************************************************************
// Load Function at the tip of the bar
//***************************************************************************************************************************************************
double LoadFunction ( double& Time, double& Alpha, double& P) 
{
double LoadFac ;

if (Time < 2*pi/Alpha) LoadFac = -P * sin( Alpha * Time );
else LoadFac = 0.0;

  return (LoadFac);
}

//***************************************************************************************************************************************************
// computing exact solution
//***************************************************************************************************************************************************
void HistorySolution ( int& NJ, double& Time, double& Alpha, double& P, double& E, double& Rho, double& A,   double *& U_EX, double **& XYZ) 
{

double fac; 
double x; 
double c; 

c = sqrt(E/Rho);

fac = P * c / (Alpha * E * A);

  for (int ij=0;ij<NJ;ij++) {
    x = XYZ[ij][0];
    U_EX[ij] = 0.0;

    if ( Time >= (x/c)                  ) U_EX[ij] = U_EX[ij] + fac * (1 - cos( Alpha * (Time-x/c) ));
    if ( Time >= (x/c + 2 * pi / Alpha) ) U_EX[ij] = U_EX[ij] - fac * (1 - cos( Alpha * (Time-x/c - 2 * pi / Alpha ) ));

  }

}


//***************************************************************************************************************************************************
// Domain Reduction Load  - The function creates the total load vector 
//***************************************************************************************************************************************************
void DRM_Loads_Implicit ( double & alpha1, double & alpha2, double & Time, int NDim, int NNBndry, int NNLayer, int& Wave_Type, int& Wave_Func, double & amplitude, double & c, double *& UN, double **& XYZ, int *& NoBndry_DRM, int *& NoLayer_DRM, double **& M_eb, double **& C_eb, double **& K_eb, int *&ND_e, int *&ND_b )
{

// = Local Variables ================================================================================================================================
int   i,j;       // Loop index on the nodes

double x;        // The coordinate
double u ;   // Analytical displacement
double v ;   // Analytical velocity
double a ;   // Analytical acceleration

double * U_b;    // The vector that holds the analytical acceleration at the boundary nodes. In this 1D problem, there is only one node.
double * Ud_b;   // The vector that holds the analytical velocity at the boundary nodes. In this 1D problem, there is only one node. (Not needed because we do not have damping in the system)
double * Udd_b;  // The vector that holds the analytical acceleration at the boundary nodes. In this 1D problem, there is only one node.
double * F_b;    // The vector that holds the loads for the boundary nodes.

double * U_e;    // The vector that holds the analytical displacement at the boundary nodes. In this 1D problem, there is only one node.
double * Ud_e;   // The vector that holds the analytical acceleration at the boundary nodes. In this 1D problem, there is only one node.  (Not needed because we do not have damping in the system)
double * Udd_e;  // The vector that holds the analytical acceleration at the boundary nodes. In this 1D problem, there is only one node.
double * F_e;    // The vector that holds the loads for the layer nodes.


// = Function =======================================================================================================================================
//c = sqrt(E/Rho);

// Defining the required vectors
U_b   = new double [ NNBndry * NDim ];
Ud_b  = new double [ NNBndry * NDim ];
Udd_b = new double [ NNBndry * NDim ];

F_e   = new double [ NNLayer * NDim ];

  for ( i=0; i<NNLayer * NDim; i++ ) {
    F_e [i] = 0.0;
  }


// Loop on the nodes on the DRM boundary to find out the analytical solution (In this case only one node)
  for (i=0;i<NNBndry;i++) {
    for (j=0;j<NDim;j++) {
      x  = XYZ [ NoBndry_DRM[i] ][j];  // Coordinate of the node
      u = v = a = 0.0; // Initialize the values
        // Computing the analytical solution - Comment: the one-dimensional wave-motion is identical for both SV and P waves.
        if      ( Wave_Type == 0 ) DRM_PointValues ( Wave_Func, amplitude, Time, x, c, omega, alpha1, alpha2, u, v, a); // SV wave
        else if ( Wave_Type == 1 ) DRM_PointValues ( Wave_Func, amplitude, Time, x, c, omega, alpha1, alpha2, u, v, a); // P wave

      // Filling the analytical solution vector
      U_b  [i * NNBndry * NDim + j] = u ;
      Ud_b [i * NNBndry * NDim + j] = v ;
      Udd_b[i * NNBndry * NDim + j] = a ;
    }
  }

  // Multiply the Mass, Damp, and Stiffness matrix by the vector
  for (i=0;i<NNLayer*NDim;i++) {
    for (j=0;j<NNBndry*NDim;j++) {
      F_e [i] += M_eb[i][j] * Udd_b [j] + C_eb[i][j] * Ud_b [j] + K_eb[i][j] * U_b [j]  ;
    }
  }

// Assemble the load vector 
  for (i=0; i <NNLayer*NDim; i++){
   UN[ ND_e [i] ] += F_e[i];
  }

delete U_b  ;
delete Ud_b ;
delete Udd_b;
delete F_e  ;


U_e   = new double [ NNLayer * NDim ];
Ud_e  = new double [ NNLayer * NDim ];
Udd_e = new double [ NNLayer * NDim ];

F_b   = new double [ NNBndry * NDim ];

  for ( i=0; i<NNBndry * NDim; i++ ) {
    F_b [i] = 0.0;
  }

  // Loop on the nodes on the DRM layer to find out the analytical solution (In this case only two nodes)
  for (i=0;i<NNLayer; i++) {
    for (j=0;j<NDim;j++) {
      x  = XYZ [ NoLayer_DRM[i] ][j];  // Coordinate of the node
      u = v = a = 0.0; // Initialize the values
        // Computing the analytical solution - Comment: the one-dimensional wave-motion is identical for both SV and P waves.
        if      ( Wave_Type == 0 ) DRM_PointValues ( Wave_Func, amplitude, Time, x, c, omega, alpha1, alpha2, u, v, a); // SV wave
        else if ( Wave_Type == 1 ) DRM_PointValues ( Wave_Func, amplitude, Time, x, c, omega, alpha1, alpha2, u, v, a); // P wave

      // Filling the analytical solution vector
      U_e  [i * NNBndry * NDim + j] = u ;
      Ud_e [i * NNBndry * NDim + j] = v ;
      Udd_e[i * NNBndry * NDim + j] = a ;

    }
  }

  // Multiply the Mass, Damp, and Stiffness matrix by the vector
  for (i=0;i<NNBndry*NDim;i++) {
    for (j=0;j<NNLayer*NDim;j++) {
      F_b [i] += - ( M_eb[j][i] * Udd_e [j] + C_eb[j][i] * Ud_e [j] + K_eb[j][i] * U_e [j] ) ;
    }
  }

  // Assemble the load vector 
  for (i=0; i <NNBndry*NDim; i++){
   UN[ ND_b [i]  ] += F_b[i]; 
  }

delete U_e  ;
delete Ud_e ;
delete Udd_e;
delete F_b;

}

//***************************************************************************************************************************************************
// Analytical solution in the time domain 
//***************************************************************************************************************************************************
void DRM_PointValues (int & Wave_Func, double & amplitude, double & Time, double & x, double & c, double & omega, double & alpha1, double & alpha2, double & u, double & v, double & a)
{

// = Local Variables ================================================================================================================================

// - Integer variables ------------------------------------------------------------------------------------------------------------------------------
int j;       // loop index

int TotalCycle;
int Direction ;

// - Real variables ---------------------------------------------------------------------------------------------------------------------------------
double arg1; // wave phase for the moving in the positive direction
double arg2; // wave phase for the moving in the negative direction

double LowerLimit; // Lower limit of the phase
double UpperLimit; // Upper limit of the phase

double wr;         // Central frequency of the Ricker pulse
double t_max;      // Loading duration
double fr; 

double ur;

// - Strings ----------------------------------------------------------------------------------------------------------------------------------------

// - bool -------------------------------------------------------------------------------------------------------------------------------------------

// - integer arrays ---------------------------------------------------------------------------------------------------------------------------------

// - real arrays ------------------------------------------------------------------------------------------------------------------------------------

// - data structure ---------------------------------------------------------------------------------------------------------------------------------

// ==================== Code ========================================================================================================================

// The analytical solution is u (x,t) = Ui (f_inc (t - x/c) + f_inc (t + x/c) )
// phases

  switch (Wave_Func) {
    case 0:   // sine function 

      wr    = 2.0 * pi * omega ; // characteristic central circular frequency

      arg1 = Time - x/c ; // positive direction phase - incident wave
      arg2 = Time + x/c ; // negative direction phase - reflected wave

      arg1 = arg1 * wr ;
      arg2 = arg2 * wr ;

      LowerLimit = alpha1 ;
      UpperLimit = alpha2 ;

        if ( LowerLimit <= arg1 && arg1 <= UpperLimit ) {

          u = u  +           amplitude * sin(arg1) ;
          v = v  + wr *      amplitude * cos(arg1) ;
          a = a  - wr * wr * amplitude * sin(arg1) ;

        }

        if ( LowerLimit <= arg2 && arg2 <= UpperLimit ) {

          u = u  +           amplitude * sin(arg2) ;
          v = v  + wr *      amplitude * cos(arg2) ;
          a = a  - wr * wr * amplitude * sin(arg2) ;

        }

    break;
    case 1:   // Ricker

      TotalCycle = (int)(alpha1);
      Direction  = (int)(alpha2);

      fr    = omega ;      // Central frequency of Ricker pulse
      wr    = 2.0 * pi * fr ; // characteristic central circular frequency

      //arg1 = wr * Time - wr * x/c ; // positive direction phase - incident wave
      //arg2 = wr * Time + wr * x/c ; // negative direction phase - reflected wave

      arg1 = Time - x/c ; // positive direction phase - incident wave
      arg2 = Time + x/c ; // negative direction phase - reflected wave

      //wr = omega ;           // characteristic central circular frequency

      t_max = 6.0 * sqrt( 6.0) / wr ;  // duration of loading

      LowerLimit = 0.0;
      UpperLimit = t_max ;

        if ( LowerLimit <= arg1 && arg1 <= UpperLimit ) {

          ur     = wr * arg1 - 3.0 * sqrt( 6.0 ) ;

          u  = u  + amplitude * ( ( 0.25 * ur * ur - 0.5 ) * exp( -0.25 * ur * ur ) - 13.0 * exp( -13.5 ) ) / ( 0.5 + 13.0 * exp( -13.5 ) ) ;
          v  = v  + amplitude * ( wr * ( 0.75 * ur - 0.125 * pow(ur,3.0 ) )* exp( -0.25 * ur * ur ) ) / ( 0.5 + 13.0 * exp( -13.5 ) ) ;
          a  = a  + amplitude * ( pow(wr ,2.0) * ( 0.75 - 0.75 * pow(ur,2.0) + 0.0625 * pow(ur,4.0 ) ) * exp( -0.25 * ur * ur ) ) / ( 0.5 + 13.0 * exp( -13.5 ) ) ;
        }
        else {
          for (j=2;j<=TotalCycle; j++){
            arg1 = Time - x/c -(j-1)*t_max ; 
              if ( LowerLimit <= arg1 && arg1 <= UpperLimit ) {
                ur     = wr * arg1 - 3.0 * sqrt( 6.0 ) ;

                u  = u  + (pow(Direction,(j-1))) * amplitude * ( ( 0.25 * ur * ur - 0.5 ) * exp( -0.25 * ur * ur ) - 13.0 * exp( -13.5 ) ) / ( 0.5 + 13.0 * exp( -13.5 ) ) ;
                v  = v  + (pow(Direction,(j-1))) * amplitude * ( wr * ( 0.75 * ur - 0.125 * pow(ur,3.0 ) )* exp( -0.25 * ur * ur ) ) / ( 0.5 + 13.0 * exp( -13.5 ) ) ;
                a  = a  + (pow(Direction,(j-1))) * amplitude * ( pow(wr ,2.0) * ( 0.75 - 0.75 * pow(ur,2.0) + 0.0625 * pow(ur,4.0 ) ) * exp( -0.25 * ur * ur ) ) / ( 0.5 + 13.0 * exp( -13.5 ) ) ;
              }
          }
        }

        if ( LowerLimit <= arg2 && arg2 <= UpperLimit ) {

          ur     = wr * arg2 - 3.0 * sqrt( 6.0 ) ;

          u  = u  + amplitude * ( ( 0.25 * ur * ur - 0.5 ) * exp( -0.25 * ur * ur ) - 13.0 * exp( -13.5 ) ) / ( 0.5 + 13.0 * exp( -13.5 ) ) ;
          v  = v  + amplitude * ( wr * ( 0.75 * ur - 0.125 * pow(ur,3.0 ) )* exp( -0.25 * ur * ur ) ) / ( 0.5 + 13.0 * exp( -13.5 ) ) ;
          a  = a  + amplitude * ( pow(wr ,2.0) * ( 0.75 - 0.75 * pow(ur,2.0) + 0.0625 * pow(ur,4.0 ) ) * exp( -0.25 * ur * ur ) ) / ( 0.5 + 13.0 * exp( -13.5 ) ) ;
        }
        else {
          for (j=2;j<=TotalCycle; j++){
            arg2 = Time + x/c -(j-1)*t_max;
              if ( LowerLimit <= arg2 && arg2 <= UpperLimit ) {

                ur     = wr * arg2 - 3.0 * sqrt( 6.0 ) ;

                u  = u  + (pow(Direction,(j-1))) * amplitude * ( ( 0.25 * ur * ur - 0.5 ) * exp( -0.25 * ur * ur ) - 13.0 * exp( -13.5 ) ) / ( 0.5 + 13.0 * exp( -13.5 ) ) ;
                v  = v  + (pow(Direction,(j-1))) * amplitude * ( wr * ( 0.75 * ur - 0.125 * pow(ur,3.0 ) )* exp( -0.25 * ur * ur ) ) / ( 0.5 + 13.0 * exp( -13.5 ) ) ;
                a  = a  + (pow(Direction,(j-1))) * amplitude * ( pow(wr ,2.0) * ( 0.75 - 0.75 * pow(ur,2.0) + 0.0625 * pow(ur,4.0 ) ) * exp( -0.25 * ur * ur ) ) / ( 0.5 + 13.0 * exp( -13.5 ) ) ;
              }
          }
        }

    break;
    default:
      cout << "Wave function is not defined."<< endl;
  }

}

//***************************************************************************************************************************************************
// Analytical solution in the frequency domain 
//***************************************************************************************************************************************************
void DRM_PointValues_Freq ( double & amplitude, double & x, double & c, double & omega, double & u_R, double & u_I )
{

// = Local Variables ================================================================================================================================

// - Integer variables ------------------------------------------------------------------------------------------------------------------------------
//int j;       // loop index

// - Real variables ---------------------------------------------------------------------------------------------------------------------------------

double k;    // wavenumber

// - Strings ----------------------------------------------------------------------------------------------------------------------------------------

// - bool -------------------------------------------------------------------------------------------------------------------------------------------

// - integer arrays ---------------------------------------------------------------------------------------------------------------------------------

// - real arrays ------------------------------------------------------------------------------------------------------------------------------------

// - data structure ---------------------------------------------------------------------------------------------------------------------------------

// ==================== Code ========================================================================================================================

// The analytical solution is u (x,t) = u_i (exp(i k x) + exp (-i k x)) = 2 u_i cos(kx)

k = omega / c;  // wavenumber

u_R = 2 * amplitude * cos (k * x);  // The real part of the analytical solution in the frequency domain
u_I = 0.0;                          // The imaginary part of the analytical solution in the frequency domain



}


