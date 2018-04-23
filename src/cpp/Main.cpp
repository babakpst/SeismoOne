/*
===============================================================================

********************G E T A F I X   1 D ***************************************

Developer:         Babak Poursartip
Code description:  Numerical solution of one-dimensional wave propagation in a heterogeneous half-space
Date starting:     July 1, 2016
version 1.0        July 07, 2016     Time domain
version 1.1        July 20, 2016     output for Matlab
version 1.2        July 25, 2016     Matlab output has modified
version 1.3        July 27, 2016     Visualization (Animation)
version 1.4        Sept 10, 2016     Unstructured mesh
version 2.0        Sept 22, 2016     Transfer functions (frequency domain)
version 2.0        March 30, 2018    Check
version 2.1        April 08, 2018    Some minor modifications
version 2.2        April 08, 2018    Transforming the code to an OOP

Last Update:       April 12, 2018

Comments:
Node number starts from 0
Element number starts from 0
Equation number starts from 0
  ID(1 2
     3 4 )
  ND [---x--- ---y---]

===============================================================================
*/

#include "../include/Variables.h"
#include "../include/Input.h"
#include "../include/Discretization.h"
#include "../include/GlobalMatrices_Full.h"
#include "../include/GlobalMatrices_Skyline.h"
#include "../include/Solver_Full.h"
#include "../include/Solver_Skyline.h"
#include "../include/Fre_Full.h"

int main()
{



// = Open input files ===============================================================================================================================
// Address file
ifstream Addressfile;
Addressfile.open ("Address.txt", ios::in );

// Reading simulation information
getline (Addressfile,TempS);
Addressfile >> Name ; // Input file name

getline (Addressfile,TempS);
getline (Addressfile,TempS);
Addressfile >> Directory ; // Directory of the input file

cout << endl;
cout << "----- Analysis information -----" << endl;
cout << " Name:  " << Name << endl;
cout << " Directory:  " << Directory << endl;

// Windows 
/*
Input_Dir            = Directory + "\\Input\\"   + Name + ".txt";
OutputMatlab_Dir     = Directory + "\\Output\\"  + Name + ".Matlab";
Info_Dir             = Directory + "\\Output\\"  + Name + ".inf";

FullFile_Dir         = Directory + "\\Output\\"  + Name + ".Res";
HistoryFile_Dir      = Directory + "\\Output\\"  + Name + ".His";
TransferFunction_Dir = Directory + "\\Output\\"  + Name + ".TRF";
*/

// Linux

Input_Dir            = Directory + "/input/"   + Name + ".txt";
OutputMatlab_Dir     = Directory + "/output/"  + Name + ".Input";
Info_Dir             = Directory + "/output/"  + Name + ".Input";

FullFile_Dir         = Directory + "/output/"  + Name + ".Input";
HistoryFile_Dir      = Directory + "/output/"  + Name + ".Input";
TransferFunction_Dir = Directory + "/output/"  + Name + ".Input";


// Final Directories
cout<< "  Input dir:   " << Input_Dir << endl;
cout<< "  Output dir:  " << Info_Dir  << endl;

// Input File
ifstream InputFile;
InputFile.open (Input_Dir.c_str(), ios::in );

// Output file for Matlab for visualization
ofstream OutputMatlab;
OutputMatlab.open (OutputMatlab_Dir.c_str(), ios::out );

// Information file
ofstream info;
info.open (Info_Dir.c_str(), ios::out );

// Reading basic data
// = Reading data ===================================================================================================================================
InputBasic (NEl_DRM, Solver, Dis_History, alpha1, alpha2, amplitude, omega, Wave_Type, Wave_Func, NInt, NDim, NNode, NDOF, NMat, NPM, NNodePWaveL, LoadType, Beta, Gama, DT, NStep, L, Alpha, P, A, OShFunc, InputFile);

// Allocating required 1D arrays - vectors
Length      = new double[NMat];  // Material Type of Elements
NoBndry_DRM = new int[NNBndry];  // vector that holds the node numbers on the DRM boundary- There is only one node on the DRM boundary
NoLayer_DRM = new int[NNLayer];  // vector that holds the node numbers on the DRM layer. Ther are only two nodes in the layer in a 1D wave motion.

Nodal_History = new int[Dis_History];
Loc_History = new double[Dis_History];

PMat = new double *[NMat];  // Properties of Materials
  for (int i = 0; i<NMat; i++){
    PMat[i] = new double[NPM];
  }

Element_Layer = new int[NMat];
Layer_Depth   = new int[NMat];

InputArrays( omega, NNodePWaveL, Dis_History, NMat, Length, L, Loc_History, Element_Layer, Layer_Depth, NEl, OShFunc, NJ, PMat, InputFile);


// Allocating required 1D arrays - vectors
MTel = new int [ NEl ];  // Material Type of Elements

// Allocating required 2D arrays - matrices
INod = new int*[NNode];  // node connectivity
  for(int i=0;i<NNode;i++){
    INod[i]=new int[NEl];
  }

ID   = new int *[ NJ ];  // Identifications
  for(int i=0;i<NJ;i++){
    ID[i]=new int[NDOF];
  }

XYZ  = new double *[ NJ ];  // Coordinates
  for(int i=0;i<NJ;i++){
    XYZ[i]=new double[NDim];
  }

// discretization
Discretization( NEl_DRM, NDOF, Dis_History, NDim, NMat, NJ, OShFunc, NEl,  NEqM,           L,                  Length,   MTel, INod, ID, XYZ, Element_Layer, Layer_Depth, NoBndry_DRM, NoLayer_DRM, Nodal_History, Loc_History     );

// Write information
info << "Number of elements :" << NEl << endl; 
info << "Number of nodes    :" << NJ << endl; 


// Creating the input file for visualization code in Matlab
InputMatlab ( NJ, Dis_History, NMat, NEl, DT, NStep, L, OutputMatlab, NEqM, NPM, Length, PMat, Loc_History, XYZ ) ;



// Close files 
// Address file
//Addressfile.close();
InputFile.close();
OutputMatlab.close();

// Allocating required arrays
  switch (Solver) {
    case 0:   // Full matrices
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
    case 1:   // Skyline method
      JD   = new int [NEqM] ;
      NTK  = new int [NEqM] ;

      Skyline ( NEqM, NEl, NNode, NDOF, NTK, INod, ID, JD );

      K_S  = new double [ JD[NEqM-1] ];  // Stiffness Matrix
      C_S  = new double [ JD[NEqM-1] ];  // Damping matrix
      M_S  = new double [ JD[NEqM-1] ];  // Mass matrix
    break;
    case 2: // 
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

K_eb  = new double *[ NDim * NNLayer ];  // 
  for(int i=0;i<(NDim * NNLayer);i++){
    K_eb[i]=new double[ NDim * NNBndry ];
  }

C_eb  = new double *[ NDim * NNLayer ];  // 
  for(int i=0;i<(NDim * NNLayer);i++){
    C_eb[i]=new double[ NDim * NNBndry ];
  }

M_eb  = new double *[ NDim * NNLayer ];  // 
  for(int i=0; i<(NDim * NNLayer); i++){
    M_eb[i]= new double[ NDim * NNBndry ];
  }

ND_b = new int  [ NNBndry * NDim ];
ND_e = new int  [ NNLayer * NDim ];

F  = new double [NEqM] ;

  // Filling the index for layered nodes
  for (int i=0;i<NNLayer;i++) {
    for ( int j=0;j<NDim;j++) {
      ND_e [ j * NNLayer + i ] = ID [ NoLayer_DRM [ i ] ][j];
    }
  }

  // Filling the index for boundary nodes
  for ( int i=0;i<NNBndry;i++) {
    for (int j=0;j<NDim;j++) {
      ND_b [ j * NNBndry + i ] = ID [ NoBndry_DRM[i]][j];
    }
  }

// Allocating required arrays
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
  return 0;
}


