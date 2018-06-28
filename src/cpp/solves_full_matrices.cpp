

#include "../include/solves_full_matrices.h"


/*
###################################################################################################
Purpose: This function solves the AX=B using Gaussian-Elemination method.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/28/2018 - Subroutine initiated.
V0.01: 06/28/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/
void main_ns::Solver_ns::solve_full_matrices_cls::Gaussian ( int& NEqM, double *& UN, double **& K)
{

int k, l; 
double Fac;
double sum;

// Uppertriangular
std::cout << "Upper" << std::endl;
  for (int j=0; j<NEqM-1; j++) {
    for (int i=j+1; i<NEqM-1; i++) {
      Fac = -K[i][j] / K[j][j];
          for (k=j; k<NEqM; k++) {
            K[i][k] += Fac* K[j][k];
          }
        UN[i] += Fac * UN[j];
    }  
  }

//Backwardsubstitution
std::cout << "Backward" << std::endl;
  for (int i=0; i<NEqM; i++) {
    k = NEqM-i-1;
    sum = 0.0;
      for (int j=0; j<i; j++) {
        l=NEqM-j-1;
        sum += K[k][l] * UN[l];
      }
    UN[k] = (UN[k] - sum) / K[k][k];
  }
}


/*
###################################################################################################
Purpose: This function solves the AX=B using LDLT method.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/28/2018 - Subroutine initiated.
V0.01: 06/28/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/
void main_ns::Solver_ns::solve_full_matrices_cls::LDLT(int &NEqM, double **&K)
{

  int tempI;
  double *L;

  L = new double[NEqM]; // Identifications

  for (int j = 0; j < NEqM; j++)
  {
    std::cout << j << " out of " << NEqM << std::endl;

    tempI = 5 + j;
    if (tempI > NEqM)
      tempI = NEqM;

    for (int i = j + 1; i < tempI; i++)
    {
      L[i] = K[i][j] / K[j][j];
    }
    for (int k = j + 1; k < tempI; k++)
    {
      for (int l = j + 1; l < NEqM; l++)
      {
        K[k][l] = K[k][l] - L[k] * K[j][l];
      }
    }
    for (int i = j + 1; i < tempI; i++)
    {
      K[i][j] = L[i];
    }
  }
}

/*
###################################################################################################
Purpose: This function solves the AX=B using LDLT method.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/28/2018 - Subroutine initiated.
V0.01: 06/28/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/

void main_ns::Solver_ns::solve_full_matrices_cls::Substitute(int &NEqM, double *&UN, double **&K)
{

  int k, l; // Loop index
  double temp;
  double *L;

  L = new double[NEqM]; // Identifications

  //cout << "Forward" << endl;
  for (int i = 0; i < NEqM; i++)
  {
    temp = 0.0;
    for (int j = 0; j < i; j++)
    {
      temp += K[i][j] * UN[j];
    }
    UN[i] = UN[i] - temp;
  }

  for (int i = 0; i < NEqM; i++)
  {
    UN[i] = UN[i] / K[i][i];
  }

  //cout << "Backward" << endl;
  for (int i = 0; i < NEqM; i++)
  {

    k = NEqM - i - 1;
    temp = 0.0;
    for (int j = 0; j < i; j++)
    {
      l = NEqM - j - 1;
      temp += K[l][k] * L[l];
    }
    L[k] = (UN[k] - temp);
  }

  for (int i = 0; i < NEqM; i++)
  {
    UN[i] = L[i];
  }
}

/*
###################################################################################################
Purpose: This function solves the time domain problem using the Newmark method.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	
================================= V E R S I O N ===================================================
V0.00: 06/28/2018 - Subroutine initiated.
V0.01: 06/28/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/

void main_ns::Solver_ns::solve_full_matrices_cls::
     solve_the_system_using_implicit_newmark_method(
       
                                                    double &L, int &Wave_Type, int &Wave_Func, 
                                                    int &NStep, int &NEqM, int &LoadType, 
                                                    double &Alpha, double **&M, double **&C, 
                                                    double **&K, double *&F, double **&PMat, 
                                                    double **&XYZ, ofstream &FullSol, 
                                                    ofstream &History, int *&ND_e, int *&ND_b, 
                                                    int *&Nodal_History)
{

  int ij; // Loop indices
  
  // Newmark constants
  double A0 = 1.0 / (Model->Beta * Model->DT * Model->DT);
  double A1 = Model->Gama / (Model->Beta * Model->DT);
  double A2 = 1.0 / (Model->Beta * Model->DT);
  double A3 = 1.0 / (2.0 * Model->Beta) - 1.0;
  double A4 = Model->Gama / Model->Beta - 1.0;
  double A5 = Model->DT * (Model->Gama / (2.0 * Model->Beta) - 1.0);




  double Time;                   // time
  double Total_Time;             // time

  double TE;                     // temporary variable
  double Initial_Time; // Starting time of the simulation

  double *UN;   // temporay arrays
  double *U;    // temporay arrays
  double *UD;   // temporay arrays
  double *UDD;  // temporay arrays
  double *Temp; // temporay arrays

  bool InitialTime = false;


  double E = PMat[0][0];    // Elastic Modulus of the base material required for the DRM loads
  double Rho = PMat[0][1];  // density of the base material required for the DRM loads
  double c = sqrt(E / Rho); // wave velocity of the base material required for the DRM loads
  
  
  

  UN   = new double[NEqM];
  U    = new double[NEqM];
  UD   = new double[NEqM];
  UDD  = new double[NEqM];
  Temp = new double[NEqM];



  // Initializing displacement, velocity and acceleration
  for (int i = 0; i < NEqM; i++)
  {
    UN[i] = 0.0;
    U[i] = 0.0;
    UD[i] = 0.0;
    UDD[i] = 0.0;
  }

  // Effective stiffness matrix
  std::cout << "Effective stress ..." << std::endl;
  for (int i = 0; i < NEqM; i++)
  {
    for (int j = 0; j < NEqM; j++)
    {
      K[i][j] = K[i][j] + A0 * M[i][j] + A1 * C[i][j];
    }
  }

  // Reduction the coefficient matrix ( Effective Stiffness Matrix )
  //Reduce_Full (NEqM, K, Check);
  LDLT(NEqM, K);

  Initial_Time = -L / c;

  // Solve the PDE for each step -----------------------------------------------
  for (int IStep = 0; IStep < NStep + 1; IStep++)
  {

    Total_Time = IStep * DT;
    Time = Initial_Time + IStep * DT; // Time STEP
    std::cout << "Time Step:  " << IStep << "  Time: " << Time << "  Total time:" << Total_Time << std::endl;

    // Update displacements, velocity and acceleration
    for (int i = 0; i < NEqM; i++)
    {
      Temp[i] = UD[i];
      U[i] = UN[i] - U[i];
      UD[i] = A1 * U[i] - A4 * Temp[i] - A5 * UDD[i];
      UDD[i] = A0 * U[i] - A2 * Temp[i] - A3 * UDD[i];
      U[i] = UN[i];
      UN[i] = 0.0;
    }

    //
    //// Recording the full results
    //FullSol << "Time:  " << Time << "\n";
    //HistorySolution ( NJ, Time, Alpha, P, E, Rho, A,   U_EX, XYZ) ;
    //  // full results
    //  for (ij=0;ij<NJ;ij++) {
    //    if ( ID [ij][0] < 0 ) u = 0.0;
    //    else u = U[ij];
    //    FullSol << setw(15) << XYZ[ij][0] << setw(15) << u << setw(15) << U_EX[ij] << "\n";
    //  }
    //

    // Effective force - stored in UN
    for (int i = 0; i < NEqM; i++)
    { // find the coefficient of the M matrix
      Temp[i] = A0 * U[i] + A2 * UD[i] + A3 * UDD[i];
    }
    for (int i = 0; i < NEqM; i++)
    { // multiplying the vector by the mass matrix
      TE = 0.0;
      for (int j = 0; j < NEqM; j++)
      {
        TE += M[i][j] * Temp[j];
      }
      UN[i] = TE;
    }

    for (int i = 0; i < NEqM; i++)
    {
      Temp[i] = A1 * U[i] + A4 * UD[i] + A5 * UDD[i];
    }

    for (int i = 0; i < NEqM; i++)
    {
      TE = 0.0;
      for (int j = 0; j < NEqM; j++)
      {
        TE += C[i][j] * Temp[j];
      }
      UN[i] += TE;
    }

    // Adding load at this time step
    if (LoadType == 0)
    {
      LoadFactor = LoadFunction(Time, Alpha, P); // Pressure load
      for (ij = 0; ij < NJ; ij++)
      {
        UN[ij] = UN[ij] - F[ij] * LoadFactor;
      }
    }
    else if (LoadType == 1)
    { //        DRM_Load ();
      F[0] = 0;

      DRM_Loads_Implicit(alpha1, alpha2, Time, NDim, NNBndry, NNLayer, Wave_Type, Wave_Func, amplitude, c, UN, XYZ, NoBndry_DRM, NoLayer_DRM, M_eb, C_eb, K_eb, ND_e, ND_b);
    }

    // Check whether the initial time is small enough
    if (IStep == 0)
    {
      for (int i = 0; i < NEqM; i++)
      {
        if (UN[i] != 0)
        {
          InitialTime = true;
          break;
        }
      }
      if (InitialTime == true)
      {
        std::cout << "WARNING: REDUCE THE INITIAL TIME" << std::endl;
        return;
      }
    }

    // SOLVE
    //Gauss_El_Full ( NEqM, UN, K);
    //Gaussian ( NEqM, UN, K);
    Substitute(NEqM, UN, K);

    // time history of the solution at some particular nodes
    History << setw(6) << Total_Time;
    for (int i = 0; i < Dis_History; i++)
    {
      History << setw(20) << UN[Nodal_History[i]];
    }
    History << std::endl;

    //FullSol << "Time= " << Time << endl;
    for (int i = 0; i < NEqM; i++)
    {
      FullSol << UN[i] << setw(20);
    }
    FullSol << std::endl;
  }
}
