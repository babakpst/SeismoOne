

#include "../include/ShapeFunctions_SecondOrder_cls.h"

main_ns::ShapeFunctions_ns::ShapeFunctions_SecondOrder_cls::ShapeFunctions_SecondOrder_cls(){};


/*
###################################################################################################
Purpose: This function computes the value of the shape functions of a second order 
element (3-noded) at the x1 point.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	

================================= V E R S I O N ===================================================
V0.00: 05/14/2018 - Subroutine initiated.
V0.01: 06/02/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/
void ShapeFunctions ()  //ShapeFunc_1D_3N
{

F0 = 0.5 * x1 * (x1 - 1);
F1 = (1 - x1*x1);
F2 = 0.5 * x1 * (x1 + 1);

}

/*
###################################################################################################
Purpose: This function computes the value of the derivative of shape functions of a seconde order 
element (3-noded) at the x1 point.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	

================================= V E R S I O N ===================================================
V0.00: 05/14/2018 - Subroutine initiated.
V0.01: 06/02/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/
void DifferentialOfShapeFunctions (double& x1, double& DF0, double& DF1, double& DF2)  //Dif_ShapeFunc_1D_3N
{

DF0 = x1 - 0.5;
DF1 = -2.0 * x1;
DF2 = x1 + 0.5;

}

