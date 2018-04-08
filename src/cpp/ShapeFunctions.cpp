

#include "../header/ShapeFunctions.h"

using namespace std;


//***************************************************************************************************************************************************
// computing first order shape functions
//***************************************************************************************************************************************************
void ShapeFunc_1D_2N (double& x1, double& F0, double& F1) 
{

F0 = 0.5 * (1 - x1);
F1 = 0.5 * (1 + x1);

}

//***************************************************************************************************************************************************
// computing differentials of first order shape functions
//***************************************************************************************************************************************************
void Dif_ShapeFunc_1D_2N (double& DF0, double& DF1) 
{

DF0 = -0.5;
DF1 = +0.5;

}

//***************************************************************************************************************************************************
// computing first order shape functions
//***************************************************************************************************************************************************
void ShapeFunc_1D_3N (double& x1, double& F0, double& F1, double& F2) 
{

F0 = 0.5 * x1 * (x1 - 1);
F1 = (1 - x1*x1);
F2 = 0.5 * x1 * (x1 + 1);

}

//***************************************************************************************************************************************************
// computing differentials of first order shape functions
//***************************************************************************************************************************************************
void Dif_ShapeFunc_1D_3N (double& x1, double& DF0, double& DF1, double& DF2) 
{

DF0 = x1 - 0.5;
DF1 = -2.0 * x1;
DF2 = x1 + 0.5;

}


//***************************************************************************************************************************************************
// computing Gauss points
//***************************************************************************************************************************************************
void GAUSS_Quad_POINTS ( int& NInt, Gauss& GAUSS_POINTS)
{

  switch (NInt) {
    case 1:
      GAUSS_POINTS.XInt[0] =  0.0;  // ABSCISSAE
      GAUSS_POINTS.WInt[0] = +2.0; // WEIGHTS
    break;
    case 2:
      GAUSS_POINTS.XInt[0] = - sqrt(1.0/3.0) ; // ABSCISSAE
      GAUSS_POINTS.XInt[1] = + sqrt(1.0/3.0) ; // ABSCISSAE
      GAUSS_POINTS.WInt[0] = +1.0 ; // WEIGHTS
      GAUSS_POINTS.WInt[1] = +1.0 ; // WEIGHTS
    break;
    case 3:
      GAUSS_POINTS.XInt[0] = -sqrt(3.0/5.0);  // ABSCISSAE
      GAUSS_POINTS.XInt[1] = 0.0;             // ABSCISSAE
      GAUSS_POINTS.XInt[2] = +sqrt(3.0/5.0);  // ABSCISSAE
      GAUSS_POINTS.WInt[0] = 5.0/9.0;         // WEIGHTS
      GAUSS_POINTS.WInt[1] = 8.0/9.0;         // WEIGHTS
      GAUSS_POINTS.WInt[2] = 5.0/9.0;         // WEIGHTS
    break;
    default:
      cout << "This NInt is not available in code."<< endl;
  }

}

