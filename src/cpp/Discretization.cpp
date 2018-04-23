
// include files
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
using namespace std;

// include voids
//#include "Input.h"
#include "../include/Discretization.h"

//***************************************************************************************************************************************************
// Discretization
//***************************************************************************************************************************************************
void Discretization( int & NEl_DRM, int& NDOF, int & Dis_History, int& NDim, int& NMat, int& NJ, int& OShFunc, int& NEl,  int& NEqM,           double& L,                  double*& Length,   int*& MTel, int**& INod, int**& ID, double**& XYZ, int *& Element_Layer, int *& Layer_Depth, int*& NoBndry_DRM, int*& NoLayer_DRM, int*& Nodal_History, double*& Loc_History     )
{

// = Variables ======================================================================================================================================

// - Integer variables ------------------------------------------------------------------------------------------------------------------------------
int i, j, ij;           // loop variable
int iel;          // loop variable
int Mat_index;          // material counter
int NJ_Count;           // counter on the number of joints
int NJ_Layer;           // number of joints in the layer
int NEl_Layer;          // number of elements in the layer

// - Real variables -------------------------------------------------------------------------------------------------------------------------------
double h;               // element size
double Coordinate;      // Holds the coordinate of the node
double Level;           // 
double tol; 

// - Strings --------------------------------------------------------------------------------------------------------------------------------------
string TempS;   // Temporary variable for reading strings from input files

// - bool -----------------------------------------------------------------------------------------------------------------------------------------

// - arrays ---------------------------------------------------------------------------------------------------------------------------------------

// - data structure -------------------------------------------------------------------------------------------------------------------------------

// ==================== Code ========================================================================================================================

std::cout << "Discretization ..." << endl;

// create the coordinates of nodes
std::cout << "Coordinates ..." << endl;
NJ_Count = 0;
XYZ[0][0] = -L;
  for (i = 0; i < NMat; i++) {
    NEl_Layer = Element_Layer[i];
    NJ_Layer = OShFunc * NEl_Layer + 1;
    h = 0.50 * Layer_Depth[i] / NEl_Layer;
    h = (floor(h*10000))/10000;
      if (i == 0)
        Level = -L;
      else
        Level = Length[i-1];

      for (ij = 1; ij < NJ_Layer-2; ij++) {   // mind the error in saving real numbers
        NJ_Count = NJ_Count + 1;
        XYZ[NJ_Count][0] = Level + ij * h;      
      }
    NJ_Count = NJ_Count + 1;
    XYZ[NJ_Count][0] = Length[i] - 0.5 * (Length[i] - XYZ[NJ_Count-1][0] );
    NJ_Count = NJ_Count + 1;
    XYZ[NJ_Count][0] = Length[i];
  }

  if ((NJ-1) != NJ_Count) {
    std::cout << " Fatal Error: refer to the discretization function." << endl;
    return;
  }
  
// connectivity - Note: Node numbers start from zero
std::cout << "Connectivities ..." << endl;
  for (iel=0;iel<NEl;iel++) 
    {
      if ( OShFunc == 1 )
      {
        INod[0][iel] = iel;
        INod[1][iel] = iel+1;
      }
      else if ( OShFunc == 2 )
      {
        INod[0][iel] = iel*2 + 0;
        INod[1][iel] = iel*2 + 1;
        INod[2][iel] = iel*2 + 2;
      }
  }

// constraints
std::cout << "Equation numbers ..." << endl;
  for (i=0;i<NJ;i++) {
    for (j=0;j<NDim;j++) {
      ID[i][j] = 0;
    }
  }

// computing equation number
NEqM = -1;
  for (i=0;i<NJ;i++) {
    for (j=0;j<NDOF;j++) {
      NEqM = NEqM + 1;
      ID[i][j] = NEqM;  // We already know that no node has been fixed. Equation number then starts from zero. 
    }
  }
NEqM = NEqM + 1 ;


// material property of element
std::cout << "Material properties of elements..." << endl;
Mat_index = 0;
  for (iel=0;iel<NEl;iel++) {
    Coordinate = XYZ[ INod[2][iel] ][0];
      if (Coordinate <= Length[Mat_index])
        MTel[iel] = Mat_index;
      else
      {
        Mat_index++;
        MTel[iel] = Mat_index;
      }
  }

// Finding the DRM node numbers

//NoBndry_DRM[0] = (int)(DRM_Loc /h );
//NoLayer_DRM[0] = NoBndry_DRM [0] - 2;
//NoLayer_DRM[1] = NoBndry_DRM [0] - 1;


NoBndry_DRM[0] = 2 * NEl_DRM;
NoLayer_DRM[0] = NoBndry_DRM [0] - 2;
NoLayer_DRM[1] = NoBndry_DRM [0] - 1;

tol = 0.0001;

// Finding out the node numbers for recording the history of displacement
  for (i=0;i<Dis_History;i++) {
    for (ij=0; ij<NJ; ij++){
      if ( (Loc_History[i]-tol < XYZ[ij][0]) && (XYZ[ij][0] < Loc_History[i]+tol )) {
        Nodal_History[i] = ij;
        break;
      }
      else if ( Loc_History[i] < XYZ[ij][0] ) {
        Nodal_History[i] = ij;
        break;
      }
    }

  }

cout<< " End discretization "<< endl;

}
