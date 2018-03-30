// include files

/*
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
*/

#include "Variables.h"

#ifndef SHAPEFUNCTIONS
#define SHAPEFUNCTIONS

void ShapeFunc_1D_2N (double& x1, double& F0, double& F1);
void Dif_ShapeFunc_1D_2N ( double& DF0, double& DF1);
void ShapeFunc_1D_3N (double& x1, double& F0, double& F1, double& F2);
void Dif_ShapeFunc_1D_3N (double& x1, double& DF0, double& DF1, double& DF2);
void GAUSS_Quad_POINTS ( int& NInt, Gauss& GAUSS_POINTS);

#endif
