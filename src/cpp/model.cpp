
// include files
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <iomanip>
using namespace std;


#include "../include/Input.h"


// Constructor
main_ns::input_ns::input_cls::input_cls(){
std::cout << " --------------- SeismoOne ---------------" << endl;
std::cout << endl;
};


/*
##################################################################################################

Purpose: This function reads the data file name and directories

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin
================================= V E R S I O N ===================================================
V0.00: 05/14/2018 - Subroutine initiated.
V0.01: 05/14/2018 - Initiated: Compiled without error for the first time.

##################################################################################################
*/
void main_ns::input_ns::input_cls::Address_fun(){

// Open address files 
std::cout << " opening the address file ..." << endl;
ifstream Addressfile;
Addressfile.open ("Address.txt", ios::in );

// Reading simulation information
std::cout << " Reading file information ..." << endl;
getline (Addressfile,TempS);
Addressfile >> Name ; // Input file name

getline (Addressfile,TempS);
getline (Addressfile,TempS);
Addressfile >> Directory ; // Directory of the input file

std::cout << endl;
std::cout << " Analysis information:" << endl;
std::cout << " Name:  " << Name << endl;
std::cout << " Directory:  " << Directory << endl;

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
OutputMatlab_Dir     = Directory + "/output/"  + Name + ".Matlab";
Info_Dir             = Directory + "/output/"  + Name + ".inf";

FullFile_Dir         = Directory + "/output/"  + Name + ".Res";
HistoryFile_Dir      = Directory + "/output/"  + Name + ".His";
TransferFunction_Dir = Directory + "/output/"  + Name + ".TRF";

// Final Directories
std::cout<< "  Input dir:   " << Input_Dir << endl;
std::cout<< "  Output dir:  " << Info_Dir  << endl;

Addressfile.close();

};



//***************************************************************************************************************************************************
// reading basic data from the input file
//***************************************************************************************************************************************************
void main_ns::input_ns::input_cls::InputBasic(int &  NEl_DRM, int & Solver, int & Dis_History, double & alpha1, double & alpha2, double & amplitude, double & omega, int& Wave_Type, int& Wave_Func, int& NInt, int& NDim, int& NNode, int& NDOF, int& NMat, int& NPM, int& NNodePWaveL, int& LoadType, double& Beta, double& Gama, double& DT, int& NStep, double& L, double& Alpha, double& P, double& A, int& OShFunc, ifstream& InputFile)
{
// = Variables ======================================================================================================================================

// - Integer variables ----------------------------------------------------------------------------------------------------------------------------

// - Real variables -------------------------------------------------------------------------------------------------------------------------------
double Total_Time; // Total simulation time

// - Strings --------------------------------------------------------------------------------------------------------------------------------------
string TempS;           // Temporary variable for reading strings from input files

// ==================== Code ======================================================================================================================

cout << "Input Basic ..." << endl;

getline(InputFile,TempS);
getline(InputFile,TempS);
getline(InputFile,TempS);
getline(InputFile,TempS);

InputFile >> NMat >> NPM >> NNodePWaveL >> LoadType >> Dis_History >> Solver;

NInt =3;
NDim =1;
NDOF =1;
NNode=3;

getline(InputFile,TempS);
getline(InputFile,TempS);
getline(InputFile,TempS);


InputFile >> Beta >> Gama >> DT >> Total_Time;
NStep = (int)(Total_Time / DT);

getline(InputFile,TempS);
getline(InputFile,TempS);
getline(InputFile,TempS);

InputFile >> L >> Alpha >> P >> A >> OShFunc;

getline(InputFile,TempS);
getline(InputFile,TempS);
getline(InputFile,TempS);
InputFile >> Wave_Type >> Wave_Func >> alpha1 >> alpha2 >> amplitude >> omega >> NEl_DRM;

//cout << "Wave_Type " << Wave_Type <<" Wave_Func " <<  Wave_Func <<" alpha1 " <<  alpha1 << " alpha2 "<<alpha2 <<" amplitude "<<  amplitude << " omega " << omega<< endl;
//cin.ignore();

}

//***************************************************************************************************************************************************
// Input function for arrays
//***************************************************************************************************************************************************
void main_ns::input_ns::input_cls::InputArrays(double & omega, int & NNodePWaveL, int & Dis_History, int& NMat, double*& Length, double& L, double*& Loc_History, int*& Element_Layer, int*& Layer_Depth, int & NEl, int& OShFunc, int & NJ,  double **& PMat, ifstream& InputFile)
{

// = Variables ======================================================================================================================================

// - Integer variables ------------------------------------------------------------------------------------------------------------------------------
int imat, i;            // loop variable
int NEl_Layer;          // Number of elements in each layer

// - Real variables -------------------------------------------------------------------------------------------------------------------------------
double Layer;           // Layer size
double m1, m2;          // material properties
double Layer_check;     // Layer size
double G;               // Module (Shear or Elastic)
double rho;             // Density
double c;               // Wave velocity
double wavelength;      // wavelength of the wave
double h;               // element size

// - Strings --------------------------------------------------------------------------------------------------------------------------------------
string TempS;   // Temporary variable for reading strings from input files

// - bool -----------------------------------------------------------------------------------------------------------------------------------------

// - arrays ---------------------------------------------------------------------------------------------------------------------------------------

// - data structure -------------------------------------------------------------------------------------------------------------------------------
  
// ==================== Code ========================================================================================================================

std::cout << "Reading input files (arrays) ..." << endl;

getline(InputFile, TempS);
getline(InputFile, TempS);


std::cout << "Material properties ..." << endl;
  for (i = 0; i<NMat; i++) {
    InputFile >> imat >> m1 >> m2;
    PMat[imat - 1][0] = m1;
    PMat[imat - 1][1] = m2;
  }

std::cout << "Layeres ..." << endl;
getline(InputFile, TempS);
getline(InputFile, TempS);
Layer_check = 0.0;
  for (i = 0; i<NMat; i++) {
    InputFile >> Layer;
    Layer_check += Layer;
    Layer_Depth[i] = Layer;
    Length[i] = -L + Layer_check;
  }

  if (Layer_check != L) {  // To see if the length of layers match up with the total length
    std::cout << "Warning: Total sum of layers is not equal to the total length:  " << L << " /= " << Layer_check << endl;
    return;
  }

// Reading the location of nodes for recording the history of displacement
std::cout << "Time history locations ..." << endl;
getline(InputFile, TempS);
getline(InputFile, TempS);
  for (i = 0; i<Dis_History; i++) {
    InputFile >> Loc_History[i];
  }

// compute discretization information
NEl = 0;
  for (i = 0; i < NMat; i++) {
    G   = PMat[i][0];
    rho = PMat[i][1];
    c = sqrt(G / rho);
    wavelength = c / omega;
    h = 2 * wavelength / NNodePWaveL;  // max allowable element size
    NEl_Layer = (int)(ceil(Layer_Depth[i]/h));  // NEl_Layer
    Element_Layer[i] = NEl_Layer;
    NEl += Element_Layer[i];
  }

NJ = OShFunc * NEl + 1;

std::cout << " Total number of elements= " << NEl << endl;
std::cout << " Total number of joints  = " << NJ << endl;

std::cout << " End InputArrays " << endl;

}

//***************************************************************************************************************************************************
// This function creates input file for Matlab
//***************************************************************************************************************************************************
void main_ns::input_ns::input_cls::InputMatlab ( int & NJ, int & Dis_History, int& NMat, int& NEl, double& DT, int& NStep, double& L, ofstream& OutputMatlab, int& NEqM, int& NPM, double*& Length, double**& PMat, double*& Loc_History, double**& XYZ)
{

OutputMatlab << NMat << endl;
OutputMatlab << NEl << endl;
OutputMatlab << Dis_History << endl;
OutputMatlab << DT << endl;
OutputMatlab << NStep << endl;
OutputMatlab << L << endl;
OutputMatlab << NEqM << endl;

  for (int i = 0; i<NMat; i++){
      for (int j = 0; j<NPM; j++){
        OutputMatlab << setw(20) << PMat[i][j] ;
      }
    OutputMatlab << endl;
  }

  for (int i = 0; i<NMat; i++){
    OutputMatlab << setw(20) << Length[i] ;
  }
OutputMatlab << endl;

  for (int i = 0; i<Dis_History; i++){
    OutputMatlab << setw(20) << Loc_History [i]  ;
  }
OutputMatlab << endl;

  for (int i = 0; i<NJ; i++){
    OutputMatlab << XYZ[i][0]  << endl;
  }

OutputMatlab << endl;

}
