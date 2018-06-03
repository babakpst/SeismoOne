

#include "../include/ShapeFunctions_FirstOrder_cls.h"

main_ns::ShapeFunctions_ns::ShapeFunctions_FirstOrder_cls::ShapeFunctions_FirstOrder_cls(){};


/*
###################################################################################################
Purpose: This function computes the value of the shape functions of a first order 
element (linear- 2-noded) at the x1 point.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	

================================= V E R S I O N ===================================================
V0.00: 05/14/2018 - Subroutine initiated.
V0.01: 06/02/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/
void ShapeFunctions (double& x1, double& F0, double& F1) //ShapeFunc_1D_2N
{

F0 = 0.5 * (1 - x1);
F1 = 0.5 * (1 + x1);

}

/*
###################################################################################################
Purpose: This function computes the value of the derivative of shape functions of a first order 
element (linear- 2-noded) at the x1 point.

Developed by: Babak Poursartip
 
The Institute for Computational Engineering and Sciences (ICES)
The University of Texas at Austin	

================================= V E R S I O N ===================================================
V0.00: 05/14/2018 - Subroutine initiated.
V0.01: 06/02/2018 - Initiated: Compiled without error for the first time.

###################################################################################################
*/
void DifferentialOfShapeFunctions (double& DF0, double& DF1)  //Dif_ShapeFunc_1D_2N
{

DF0 = -0.5;
DF1 = +0.5;

}
