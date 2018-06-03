

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
template<int NInt, int NNode>
void main_ns::ShapeFunctions_ns::ShapeFunctions_cls::Retrieving_Gauss_Points_fn (const int& NInt, const int& NNode){

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


