

#include "Solver_Full.h"
using namespace std;

/*
***************************************************************************************************************************************************
 Solver: Gauss-Elimination
***************************************************************************************************************************************************
void Gaussian ( int& NEqM, double *& UN, double **& K)
{

int i, j, k, l ;         // Loop index
double Fac;
double sum;

// Uppertriangular
cout << "Upper" << endl;
  for (j=0; j<NEqM-1; j++) {
    for (i=j+1; i<NEqM-1; i++) {
      Fac = -K[i][j] / K[j][j];
          for (k=j; k<NEqM; k++) {
            K[i][k] += Fac* K[j][k];
          }
        UN[i] += Fac * UN[j];
    }  
  }


//Backwardsubstitution
cout << "Backward" << endl;
  for (i=0; i<NEqM; i++) {
    k = NEqM-i-1;
    sum = 0.0;
      for (j=0; j<i; j++) {
        l=NEqM-j-1;
        sum += K[k][l] * UN[l];
      }
    UN[k] = (UN[k] - sum) / K[k][k];
  }

}
*/


//***************************************************************************************************************************************************
// Solver: LDLT
//***************************************************************************************************************************************************
void LDLT ( int& NEqM, double **& K)
{

int i, j, k, l ;         // Loop index
int tempI;
double * L;

L   = new double [ NEqM ];  // Identifications

  for (j=0; j<NEqM; j++){
cout << j << " out of "<< NEqM<< endl;

tempI = 5+j;
if (tempI> NEqM) tempI = NEqM;


    for (i=j+1; i<tempI; i++){
      L[i] =  K[i][j] / K[j][j];
    }
    for (k=j+1; k<tempI; k++){
      for (l=j+1; l<NEqM; l++){
        K[k][l] =  K[k][l] -L[k] * K[j][l];
      }
    }
    for (i=j+1; i<tempI; i++){
      K[i][j] = L[i] ;
    }
  }

}




//***************************************************************************************************************************************************
// Solver: Substitution
//***************************************************************************************************************************************************
void Substitute ( int& NEqM, double *& UN, double **& K)
{

int i, j, k, l ;         // Loop index
double temp;
double * L;

L   = new double [ NEqM ];  // Identifications


//cout << "Forward" << endl;
  for (i=0; i< NEqM; i++) {
    temp = 0.0;
      for (j = 0; j<i; j++){
        temp += K[i][j] * UN[j];
      }
    UN[i]  = UN[i]-temp ;
  }


  for (i=0; i<NEqM; i++){
    UN[i] = UN[i]/K[i][i];
  }

//cout << "Backward" << endl;
  for (i=0; i<NEqM; i++) {

    k = NEqM-i-1;
    temp = 0.0;
      for (j=0; j<i; j++) {
       l=NEqM-j-1;
        temp += K[l][k] * L[l];
      }
    L[k] = (UN[k] - temp) ;
  }

  for (i=0; i<NEqM; i++) {
    UN[i] = L[i] ;
  }
}

//***************************************************************************************************************************************************
// Newmark's algorithm for marching in time
//***************************************************************************************************************************************************
void Newmark_Full ( double & L, int & Wave_Type, int & Wave_Func, int &NStep, int& NEqM, int& LoadType, double& Gama, double& Beta, double& DT, double& Alpha, double **& M, double **& C, double **& K, double *& F, double **& PMat, double **& XYZ, ofstream& FullSol, ofstream& History, int *&ND_e, int *&ND_b, int *&Nodal_History )
{

int    i, ij ;                  // Loop indices

double A0, A1, A2, A3, A4, A5;  // Newmark variables
double Time;                    // time
double Total_Time;              // time
double LoadFactor;              // Load Factor at each time step
double E;                       // Elastic Modulus
double Rho;                     // density
double TE ;                     // temporary variable
double c ; 
double Initial_Time ; // Starting time of the simulation

double * UN;                    // temporay arrays
double * U;                     // temporay arrays
double * UD;                    // temporay arrays
double * UDD;                   // temporay arrays
double * Temp;                  // temporay arrays

bool InitialTime = false;

E   = PMat[0][0];
Rho = PMat[0][1];
c = sqrt (E/Rho);

UN   = new double [NEqM];
U    = new double [NEqM];
UD   = new double [NEqM];
UDD  = new double [NEqM];
Temp = new double [NEqM];

// Newmark constants
A0 = 1.0 / ( Beta * DT * DT ) ;
A1 = Gama   / ( Beta * DT ) ;
A2 = 1.0 / ( Beta * DT ) ;
A3 = 1.0 / ( 2.0 * Beta ) - 1.0 ;
A4 = Gama   / Beta - 1.0 ;
A5 = DT * ( Gama / ( 2.0 * Beta ) - 1.0 ) ;

// Initializing displacement, velocity and acceleration
  for (i=0;i<NEqM;i++) {
    UN[i]  = 0.0 ;
    U[i]   = 0.0 ;
    UD[i]  = 0.0 ;
    UDD[i] = 0.0 ;
  }


// Effective stiffness matrix 
std::cout << "Effective stress ..." << endl;
  for (i=0;i<NEqM;i++) {
    for (int j=0;j<NEqM;j++) {
      K[i][j] = K[i][j] + A0 * M[i][j] + A1 * C[i][j];
    }
  }


// Reduction the coefficient matrix ( Effective Stiffness Matrix )
//Reduce_Full (NEqM, K, Check);
LDLT ( NEqM,  K);

Initial_Time = -L / c; 

  // Solve the PDE for each step -----------------------------------------------
  for (int IStep=0;IStep<NStep+1;IStep++) {

    Total_Time = IStep * DT; 
    Time  = Initial_Time + IStep * DT ;  // Time STEP
    std::cout << "Time Step:  " << IStep << "  Time: " << Time << "  Total time:" << Total_Time << endl;

    // Update displacements, velocity and acceleration
      for (i=0; i<NEqM; i++) {
        Temp[i] = UD[i] ;
        U[i]    = UN[i] - U[i] ;
        UD[i]   = A1 * U[i] - A4 * Temp[i] - A5 * UDD[i] ;
        UDD[i]  = A0 * U[i] - A2 * Temp[i] - A3 * UDD[i] ;
        U[i]    = UN[i] ;
        UN[i]   = 0.0 ; 
      }

    /*
    // Recording the full results
    FullSol << "Time:  " << Time << "\n";
    HistorySolution ( NJ, Time, Alpha, P, E, Rho, A,   U_EX, XYZ) ;
      // full results
      for (ij=0;ij<NJ;ij++) {
        if ( ID [ij][0] < 0 ) u = 0.0;
        else u = U[ij];
        FullSol << setw(15) << XYZ[ij][0] << setw(15) << u << setw(15) << U_EX[ij] << "\n";
      }
    */

      // Effective force - stored in UN
      for (i=0;i<NEqM; i++) { // find the coefficient of the M matrix
        Temp[i] = A0 * U[i] + A2 * UD[i] + A3 * UDD[i] ;
      }
      for (i=0;i<NEqM;i++) { // multiplying the vector by the mass matrix
        TE = 0.0;
          for (int j=0;j<NEqM;j++) {
            TE += M[i][j] * Temp[j];
          }
        UN[i] = TE;
      }

      for (i=0;i<NEqM; i++) {
        Temp[i] = A1 * U[i] + A4 * UD[i] + A5 * UDD[i] ;
      }

      for (i=0;i<NEqM;i++) {
        TE = 0.0;
          for (int j=0;j<NEqM;j++) {
            TE  += C[i][j] * Temp[j] ;
          }
        UN[i] +=  TE;
      }

    // Adding load at this time step
      if ( LoadType == 0 ) {
        LoadFactor = LoadFunction ( Time, Alpha, P);    // Pressure load
        for (ij=0; ij<NJ; ij++) {
          UN[ij] = UN[ij] - F[ij] * LoadFactor;
        }
      }
      else if ( LoadType == 1 ) { //        DRM_Load ();
        F[0] = 0;

        DRM_Loads_Implicit ( alpha1, alpha2, Time, NDim, NNBndry, NNLayer,              Wave_Type, Wave_Func, amplitude, c, UN, XYZ, NoBndry_DRM, NoLayer_DRM, M_eb, C_eb, K_eb, ND_e, ND_b );
      }

    // Check whether the initial time is small enough
      if ( IStep == 0) {
        for (i=0; i<NEqM; i++ ) {
          if (UN[i]!= 0) {
            InitialTime = true;
            break;
          }
        }
        if (InitialTime == true) {
          cout << "WARNING: REDUCE THE INITIAL TIME" << endl;
          return;
        }
      }

    // SOLVE
    //Gauss_El_Full ( NEqM, UN, K);
    //Gaussian ( NEqM, UN, K);
    Substitute ( NEqM, UN, K);


    // time history of the solution at some particular nodes
    History << setw(6) << Total_Time;
      for (i=0;i<Dis_History;i++) {
        History << setw(20) << UN[Nodal_History[i]]  ;
      }
    History  << endl;


    //FullSol << "Time= " << Time << endl;
      for (i = 0; i<NEqM; i++){
        FullSol << UN[i]<< setw(20);
      }
    FullSol << endl;

  }

}


