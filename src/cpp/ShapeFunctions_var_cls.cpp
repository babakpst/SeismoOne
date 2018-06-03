

#include "../include/ShapeFunctions_cls.h"

void main_ns::ShapeFunctions_ns::ShapeFunctions_cls::ShapeFunctions_cls(){};

// 
/*
###################################################################################################
Purpose: This function retrieves Gauss points for the numerical integration.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	

================================= V E R S I O N ===================================================
V0.00: 05/14/2018 - Subroutine initiated.
V0.01: 06/02/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/
void main_ns::ShapeFunctions_ns::ShapeFunctions_cls::Retrieving_Gauss_Points_fn (){

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

