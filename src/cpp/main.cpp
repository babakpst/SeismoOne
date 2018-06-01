

/*
===================================================================================================

************************************** G E T A F I X   1 D  ***************************************
Purpose:  Numerical solution of one-dimensional wave propagation in a heterogeneous half-space.

Developer:         Babak Poursartip
                   Civil, Architectural, and Environmental Engineering department
                   The Institute for Computational Engineering and Sciences (ICES)
                   The University of Texas at Austin

Date starting:     July 1, 2016
version 1.0        July 07, 2016     Time domain
version 1.1        July 20, 2016     output for Matlab
version 1.2        July 25, 2016     Matlab output has modified
version 1.3        July 27, 2016     Visualization (Animation)
version 1.4        Sept 10, 2016     Unstructured mesh
version 2.0        Sept 22, 2016     Transfer functions (frequency domain)
version 2.0        March 30, 2018    Check
version 2.1        April 08, 2018    Some minor modifications
version 2.2        May 30, 2018      Transforming the code to an OOP

Last Update:       May 20, 2018

Comments:
Node number starts from 0
Element number starts from 0
Equation number starts from 0
  ID(1 2
     3 4 )
  ND [---x--- ---y---]

===================================================================================================
*/

#include "../include/Address_cls.h"
#include "../include/Model_cls.h"
#include "../include/Discretization_cls.h"
#include "../include/Visualization_cls.h"
#include "../include/solver_cls.h"
#include "../include/Solver_Full_cls.h"

//#include "../include/Variables.h"
//#include "../include/Discretization.h"
//#include "../include/GlobalMatrices_Full.h"
//#include "../include/GlobalMatrices_Skyline.h"
//#include "../include/Solver_Full.h"
//#include "../include/Solver_Skyline.h"
//#include "../include/Fre_Full.h"

int main()
{

// = input class ==================================================================================
main_ns::address_ns::address_cls input;
input.address_fn();

// = model data ===================================================================================
main_ns::model_ns::model_cls model(&input);
model.InputBasic(); // Reading basic data
model.InputArrays(); // Reading arrays

// = discretization ===============================================================================
main_ns::discretization_ns::discretization_cls discretized_model(&model);
discretized_model.Discretization();

// = discretization ===============================================================================
main_ns::visualization_ns::visualization_cls Visual(&input, &model, &discretized_model);
Visual.MatlabOutput_fn();

// = solver =======================================================================================
main_ns::solver_ns::solver_cls* Solver;

Solver = new main_ns::solver_full_ns::solver_full_cls(&discretized_model, &model);


/*
// Allocating required arrays
  switch (Solver) {
    case 0:   // Full matrices
      Solver = new solver_full_cls(&discretized_model, &model)
    break;
    case 1:   // Skyline method
      JD   = new int [NEqM] ;
      NTK  = new int [NEqM] ;

      Skyline ( NEqM, NEl, NNode, NDOF, NTK, INod, ID, JD );

      K_S  = new double [ JD[NEqM-1] ];  // Stiffness Matrix
      C_S  = new double [ JD[NEqM-1] ];  // Damping matrix
      M_S  = new double [ JD[NEqM-1] ];  // Mass matrix
    break;
    case 2: //   Transfer functions in the frequency domain
      K  = new double *[ NEqM ];  // Stiffness Matrix
        for(int i=0;i<NEqM;i++){
          K[i]=new double[NEqM];
        }

      C  = new double *[ NEqM ];  // Damping matrix
        for(int i=0;i<NEqM;i++){
          C[i]=new double[NEqM];
        }

      M  = new double *[ NEqM ];  // Mass matrix
        for(int i=0;i<NEqM;i++){
          M[i]=new double[NEqM];
        }      
    break;
    default:
      cout << "Solver type is not available. Solver should be either 0 for full matrices or 1 for skyline method"<< endl;
  }


// Solving
  switch (Solver) {
    case 0:   // Time domain anaylsis using full matrices
      {

      // Open output files for the time domain simulations 
      ofstream FullSol;
      FullSol.open (FullFile_Dir.c_str(), ios::out );

      ofstream History;
      History.open (HistoryFile_Dir.c_str(), ios::out );

        for (int i=0; i<NEqM; i++) {
            for (int j=0; j<NEqM; j++) {
              M[i][j] = 0.0;
              C[i][j] = 0.0;
              K[i][j] = 0.0;
            }
          F[i]=0.0;
        }

      // - Computing Global Matrices -------------------------------------------------------------------------------------------------------------------
      GlobalMatrices_Full( NInt, NDim, NDOF, NNode, NEl,   MTel, INod, ID, XYZ, M, C, K ) ;

      // Get part of the global matrices for DRM loading
      DRM_Matrices_Full(  NNBndry, NNLayer, K, C, M, K_eb, C_eb, M_eb, ND_e, ND_b ) ;

      // - Newmark algorithm for marchin in time -------------------------------------------------------------------------------------------------------
      Newmark_Full ( L, Wave_Type, Wave_Func, NStep, NEqM, LoadType, Gama, Beta, DT, Alpha, M, C, K, F, PMat, XYZ, FullSol, History, ND_e, ND_b, Nodal_History ) ;

        for(int i=0;i<NEqM;i++){
          delete []K[i];
        }
      delete []K;

        for(int i=0;i<NEqM;i++){
          delete []C[i];
        }
      delete []C;

        for(int i=0;i<NEqM;i++){
          delete []M[i];
        }
      delete []M;

      // Close output files 
      FullSol.close();
      History.close();
      info.close();

      break;
    }
    case 1:   // Time domain analysis using skyline method
      {

      // Open output files for the time domain simulations 
      ofstream FullSol;
      FullSol.open (FullFile_Dir.c_str(), ios::out );

      ofstream History;
      History.open (HistoryFile_Dir.c_str(), ios::out );

        for (int i=0; i<JD[NEqM-1]; i++) {
          M_S[i] = 0.0; 
          C_S[i] = 0.0;
          K_S[i] = 0.0;
        }

        for (int i=0; i<NEqM; i++) {
          F[i]=0.0;
        }

      // - Computing Global Matrices -------------------------------------------------------------------------------------------------------------------
      GlobalMatrices_Skyline ( NInt, NDim, NDOF, NNode, MTel, INod, ID, XYZ, M_S, C_S, K_S, JD ) ;

      // Get part of the global matrices for DRM loading
      DRM_Matrices_Skyline( NNBndry, NNLayer, K_S, C_S, M_S, K_eb, C_eb, M_eb, ND_e, ND_b , JD)  ;

      // - Newmark algorithm for marchin in time -------------------------------------------------------------------------------------------------------
      Newmark_Skyline ( L, Wave_Type, Wave_Func, NStep, NEqM, LoadType, Gama, Beta, DT, Alpha, M_S, C_S, K_S, F, PMat, XYZ, FullSol, History, ND_e, ND_b, Nodal_History, JD, NTK, info );

      delete K_S;
      delete C_S;
      delete M_S;

      // Close output files 
      FullSol.close();
      History.close();
      info.close();

      break;
    }
    case 2:   // Computing transfer functions in the frequency domain anaylsis using full matrices
      {
      // Open output files for the transfer function output
      ofstream TransferFunc;
      TransferFunc.open (TransferFunction_Dir.c_str(), ios::out );

        for (int i=0; i<NEqM; i++) {
            for (int j=0; j<NEqM; j++) {
              M[i][j] = 0.0;
              C[i][j] = 0.0;
              K[i][j] = 0.0;
            }
          F[i]=0.0;
        }

      // - Computing Global Matrices -------------------------------------------------------------------------------------------------------------------
      GlobalMatrices_Full( NInt, NDim, NDOF, NNode, NEl,   MTel, INod, ID, XYZ, M, C, K ) ;

      // Get part of the global matrices for DRM loading
      DRM_Matrices_Full(  NNBndry, NNLayer, K, C, M, K_eb, C_eb, M_eb, ND_e, ND_b ) ;

      // - Computing the transfer functions in the frequency domain ------------------------------------------------------------------------------------
      Transfer_Full ( alpha1, alpha2, Wave_Type, NEqM, M, C, K, PMat, XYZ, ND_e, ND_b, TransferFunc );
        for(int i=0;i<NEqM;i++){
          delete []K[i];
        }
      delete []K;

        for(int i=0;i<NEqM;i++){
          delete []C[i];
        }
      delete []C;

        for(int i=0;i<NEqM;i++){
          delete []M[i];
        }
      delete []M;

      // Close output files 
      TransferFunc.close();
      info.close();

      break;
    }
    default:
      cout << "Solver type is not available. Solver should be either 0 for full matrices or 1 for skyline method"<< endl;
  }


// Information file
ofstream info;
info.open (Info_Dir.c_str(), ios::out );

// Write information
info << "Number of elements :" << NEl << endl; 
info << "Number of nodes    :" << NJ << endl; 



// = DeAllocating Arrays ============================================================================================================================
// DeAllocating 1D arrays
delete MTel;

// DeAllocating 2D arrays
  for(int i=0; i<NNode; i++){
   delete []INod[i];
  }
delete []INod;

  for (int i=0;i<NJ;i++){
   delete []ID[i];
  }
delete []ID;

  for(int i=0;i<NJ;i++){
   delete []XYZ[i];
  }
delete []XYZ;

// terminating the code
  cout << "End of the code \n";
  cout << "Press 'Enter' to end \n";
//  getline (cin , TempS) ;

*/

return 0;
}


